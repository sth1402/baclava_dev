#include <Rcpp.h>
#include "EPGroup_State.h"
#include "NormalizationAPI.h"
#include "Model.h"
#include "StateParameter.h"

using namespace Rcpp;

namespace {

inline double logsumexp(double a, double b) {
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  double m = std::max(a, b);
  return m + std::log1p(std::exp(-std::abs(a - b)));
}

}

void EPGroup_State::initialize_cache(const Model* model) {
  
  this->marginalize = 
    (model->I_marginalize_PRE && this->endpoint_type == Pkg::EP::PRECLINICAL) ||
    (model->I_marginalize_CEN && this->endpoint_type == Pkg::EP::CENSORED) ||
    (model->I_marginalize_CLI && this->endpoint_type == Pkg::EP::CLINICAL);
  
  constexpr uint8_t ALL_COMPONENTS = 
    static_cast<uint8_t>(Pkg::Components::H)    |
    static_cast<uint8_t>(Pkg::Components::P)    |
    static_cast<uint8_t>(Pkg::Components::BETA) |
    static_cast<uint8_t>(Pkg::Components::PSI)  |
    static_cast<uint8_t>(Pkg::Components::Z)    |
    static_cast<uint8_t>(Pkg::Components::TAU);
  
  this->logH = this->H_term(this->tau, model, ALL_COMPONENTS);
  this->logP = this->P_term(this->tau, this->Z, model, ALL_COMPONENTS);
  this->logI = this->I_term(this->tau, this->Z, model, ALL_COMPONENTS);
  this->logB = this->beta_term(this->tau, model, ALL_COMPONENTS);
  this->logCP = this->CP_term(model, ALL_COMPONENTS);
}

void EPGroup_State::initialize_cache(const Model* model, 
                                     const Rcpp::NumericVector& tau_, 
                                     const Rcpp::IntegerVector& Z_) {
  this->tau = tau_;
  this->Z = Z_;
  initialize_cache(model);
}

NumericVector EPGroup_State::calc_psi(const NumericVector& tau_eval, 
                                      const Model* model) const {
  
  NumericVector prob = model->psi[this->ipsi].parameterModel->par_from_tau_vec(tau_eval);
  prob = pmin(pmax(prob, 1e-12), 1.0 - 1e-12);
  return prob;
}

double EPGroup_State::calc_psi(const double tau_eval, 
                               const Model* model) const {
  
  double prob = model->psi[this->ipsi].parameterModel->par_from_tau(tau_eval);
  prob = std::min(std::max(prob, 1e-12), 1.0 - 1e-12);
  return prob;
}



void EPGroup_State::H_term(NumericVector& result,
                           const NumericVector& tau_eval, 
                           const Model* model, uint8_t impact, 
                           bool fallback_to_cache) const {

  /* if the parameter update is not an H component parameter or tau, skip; */
  if ((impact & static_cast<uint8_t>(Pkg::Components::H)) ||
      (impact & static_cast<uint8_t>(Pkg::Components::TAU))) {
  } else {
    if (fallback_to_cache) result = this->logH;
    return;
  }

  // if under observation at time of transition out of H, probability of 
  //   transitioning at tau_eval - t0
  // if no longer under observation at time of transition out of H, probability of 
  //   still being healthy at the time of exit (T - t0)
  double t0 = model->t0;
  const auto& H0 = model->H[0];
  
  for (int i = 0; i < this->n; ++i) {
    if (tau_eval[i] < this->endpoint_times[i]) {
      result[i] =  H0.dfunc(tau_eval[i] - t0, 0.0, /*log*/true);
    } else {
      result[i] = H0.pfunc(this->endpoint_times[i] - t0, 0.0, /*lower*/false, /*log*/true);
    }
  }

  return;
}

void EPGroup_State::P_term(NumericVector& result,
                           const NumericVector& tau_eval,
                           const IntegerVector& indolent,
                           const Model* model, uint8_t impact, bool fallback_to_cache) const {
  
  /* if the parameter update is not a P component parameter, latent Z, latent tau,
     or PSI under a marginalized model, skip; */
  if ((impact & static_cast<uint8_t>(Pkg::Components::P)) ||
      (impact & static_cast<uint8_t>(Pkg::Components::Z)) ||
      (impact & static_cast<uint8_t>(Pkg::Components::TAU)) ||
      (this->marginalize && (impact & static_cast<uint8_t>(Pkg::Components::PSI)))) {
  } else { 
    if (fallback_to_cache) result = this->logP;
    return;
  }

  result.fill(0.0);

  if (this->endpoint_type == Pkg::EP::CLINICAL) {
    
    NumericVector tau_p_hat = this->endpoint_times - tau_eval;
    result = model->P[this->iP].dfunc(tau_p_hat, tau_eval, /*log*/true);
    
    if (this->marginalize) {
      NumericVector psi = calc_psi(tau_eval, model);
      result = result + Rcpp::log1p(-psi);
    }
    result[tau_p_hat < 0.0] = R_NegInf;
    
  } else {
    double logSprog;
    const auto& Pip = model->P[this->iP];
    
    for (int i = 0; i < this->n; ++i) {
      if (tau_eval[i] > this->endpoint_times[i]) continue;
      
      double tau_p_hat = this->endpoint_times[i] - tau_eval[i];
      logSprog = Pip.pfunc(tau_p_hat, tau_eval[i], /*tail*/false, /*log*/true);
      
      if (this->marginalize) {
        double psi = calc_psi(tau_eval[i], model);
        double a = std::log(psi);
        double b = std::log1p(-psi) + logSprog;
        result[i] = logsumexp(a, b);
      } else {
        if (indolent[i] == 0) result[i] = logSprog;
      }
    }
  }
  
  return;
}

void EPGroup_State::I_term(NumericVector& result,
                           const NumericVector& tau_eval,
                           const IntegerVector& indolent, 
                           const Model* model, uint8_t impact, bool fallback_to_cache) const {
  
  /* if marginalized in Z model, skip; */
  if (this->marginalize) return;
  
  /* if the parameter update is not a PSI parameter, latent Z, or latent TAU
     with with a tau dependent psi model; */
  if ((impact & static_cast<uint8_t>(Pkg::Components::PSI)) ||
      (impact & static_cast<uint8_t>(Pkg::Components::Z)) ||
      (impact & static_cast<uint8_t>(Pkg::Components::TAU) &
       model->psi[this->ipsi].parameterModel->typeid_ != Pkg::Model::CONSTANT)){
  } else { 
    if (fallback_to_cache) result = this->logI;
    return;
  }

  result.fill(0.0);
  
  NumericVector prob = calc_psi(tau_eval, model);
  for (int i = 0; i < this->n; ++i) {
    result[i] = indolent[i] == 1 ? log(prob[i]) : log1p(-prob[i]);
  }

  return;
}

void EPGroup_State::beta_term(NumericVector& result,
                              const NumericVector& tau_eval, 
                              const Model* model, 
                              uint8_t impact, bool fallback_to_cache) const {
  
  /* if the parameter update is not a BETA or latent tau, skip; */
  if ((impact & static_cast<uint8_t>(Pkg::Components::BETA)) ||
      (impact & static_cast<uint8_t>(Pkg::Components::TAU))) {
  } else { 
    if (fallback_to_cache) result = this->logB;
    return;
  }

  result.fill(0.0);
  
  if (this->age_screen.size() == 0) return;

  NumericVector betas(model->beta.size());
  for (int i = 0; i < model->beta.size(); ++i) {
    betas[i] = model->beta[i].parameterModel->parameters[0]->value;
    betas[i] = std::max(0.0, std::min(betas[i], 1.0));
  }
  NumericVector lbetas = log(betas);
  NumericVector l1mbetas = log1p(-betas);
  
  if (this->endpoint_type == Pkg::EP::PRECLINICAL) {
    for (int i = 0; i < this->n; ++i) {
      // for each participant, determine if screens missed/caught a preclinical case
      const auto& as = this->age_screen[i];
      const auto& st = this->screen_types[i];
      int m = as.size();
      for (int j = 0; j < m-1; ++j) if(as[j] > tau_eval[i]) result[i] += l1mbetas[st[j]];
      // if preclinical, last screen caught disease
      result[i] += lbetas[st[m-1]];
    }
  } else {
    for (int i = 0; i < this->n; ++i) {
      // for each participant, determine if screens missed/caught a preclinical case
      const auto& as = this->age_screen[i];
      const auto& st = this->screen_types[i];
      int m = as.size();
      for (int j = 0; j < m; ++j) if(as[j] > tau_eval[i]) result[i] += l1mbetas[st[j]];
    }
  }

  return;
}

void EPGroup_State::CP_term(NumericVector& result, const Model* model, 
                            uint8_t impact, bool fallback_to_cache) const {

  /* Normalization is marginalized over latent tau and Z; update only when model
   * parameters are modified */
  if ((impact & static_cast<uint8_t>(Pkg::Components::PSI)) ||
      (impact & static_cast<uint8_t>(Pkg::Components::P)) ||
      (impact & static_cast<uint8_t>(Pkg::Components::H))){
  } else { 
    if (fallback_to_cache) result = this->logCP;
    return;
  }

  result = -Normalization_compute_cp_log(
    this->age_entry,
    model->t0, model->H[0], model->P[this->iP], model->psi[this->ipsi]);
  
  if (is_true(any(!is_finite(result))))
    throw NegativeRateError("normalization integral is NaN");

  return;
}

NumericVector EPGroup_State::pfunc_ab_optimized(const std::vector<double>& dist_scales, 
                                               Pkg::Dist dist_type, 
                                               const double shape_val, const double t0, 
                                               const NumericVector& wgt) const {
  
  if (dist_type != Pkg::Dist::WEIBULL && 
      dist_type != Pkg::Dist::GAMMA) {
    throw std::logic_error("pfunc_ab_optimized is only defined for continuous distributions (weibull/gamma).");
  }
  
  NumericVector result_probs(this->n_intvl);
  
  if (!(shape_val > 0.0) || !R_finite(shape_val)) {
    result_probs.fill(NA_REAL);
    return result_probs;
  }
  
  double (*p_dist_func)(double, double, double, int, int);
  
  if (dist_type == Pkg::Dist::WEIBULL) {
    p_dist_func = R::pweibull; 
  } else if (dist_type == Pkg::Dist::GAMMA) {
    p_dist_func = R::pgamma;
  } else {
    stop("Unknown distribution type.");
  }
  
  int cnt = 0; // index for the final flat result vector
  
  for (int i = 0; i < this->n; ++i) {
    const auto& ep = this->ep_values[i];
    const double et = this->endpoint_times[i];
    const double scale = dist_scales[i]; 
    
    // pre-calculate all F(t) values for this inner vector
    std::vector<double> F_values(ep.size());
    if (this->endpoint_type != Pkg::EP::CLINICAL) {
      for(int k = 0; k < ep.size(); ++k) {
        F_values[k] = p_dist_func(ep[k] - t0, shape_val, scale, /*lower=*/true, /*log=*/true);
      }
    } else {
      for(int k = 0; k < ep.size(); ++k) {
        F_values[k] = p_dist_func(et - ep[k], shape_val, scale, /*lower=*/true, /*log=*/true);
      }
    }
    
    double logFA, logFB;
    double sum_prob = 0.0;
    int cnt_start = cnt;
    int cnt_end = cnt_start;
    
    for (int j = 0; j < ep.size() - 1; ++j) {
      if (this->endpoint_type != Pkg::EP::CLINICAL) {
        logFA = F_values[j];
        logFB = F_values[j+1];
      } else {
        logFA = F_values[j+1];
        logFB = F_values[j];
      }
      
      bool bad = (!std::isfinite(logFA) && !std::isfinite(logFB)) || (logFB <= logFA);
      double diff = bad ? NA_REAL : (logFA - logFB);
      
      double logInterval = logFB + std::log1p(-std::exp(diff));
      
      if (bad || !std::isfinite(logInterval)) logInterval = R_NegInf;
      
      result_probs[cnt] = std::exp(logInterval) * wgt[cnt];
      sum_prob += result_probs[cnt];
      cnt_end = cnt;
      cnt++;
    }
    
    if (!std::isfinite(sum_prob)) throw std::invalid_argument("tau probability is not finite");
    if (sum_prob <= 0.0) throw std::invalid_argument("tau probability is not positive");
    for (int j = cnt_start; j <= cnt_end; ++j) result_probs[j] /= sum_prob;
    
  }
  
  return result_probs;
}