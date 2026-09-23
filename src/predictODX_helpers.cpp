// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "gauss_legendre_32.h"
#include "Model.h"
#include "NormalizationEngine.h"
#include "ParameterModel.h"
#include "StateParameter.h"
#include "Pkg_types.h"

using namespace Rcpp;
using gl32::gauss_legendre_32;


namespace {
inline double clamp01(double p) {
  if (p <= 1e-12) return 1e-12;
  if (p >= 1.0 - 1e-12) return 1.0 - 1e-12;
  return p;
}
  
  // piecewise-constant hazard survival S_other(age_min -> age)
  // (extends last rate beyond last tabulated age)
  double other_surv(double age_min,
                    double age,
                    const std::vector<double>& ages,
                    const std::vector<double>& rates) {
    if (!R_finite(age) || age <= age_min) return 1.0;
    
    const double am = age_min;
    const double t  = age;
    const double max_age = ages.back();
    const double t_tab = std::min(t, max_age);
    
    double hh = 0.0;
    double left = am;
    
    while (left < t_tab) {
      auto it = std::upper_bound(ages.begin(), ages.end(), left);
      int idx = std::max(0, (int)(it - ages.begin()) - 1);
      double right = t_tab;
      if (idx+1 < (int)ages.size()) right = std::min(right, ages[idx+1]);
      if (right > left) hh += rates[idx] * (right - left);
      left = right;
    }
    if (t > max_age) hh += rates.back() * (t - max_age);
    
    return std::exp(-hh);
  }
  
  double calc_psi(const double tau_eval, Model* model) {
    
    double prob = model->psi[0].parameterModel->par_from_tau(tau_eval);
    prob = std::min(std::max(prob, 1e-12), 1.0 - 1e-12);
    return prob;
    
  }
  
  double indolent_kernel(double u, Model* model) {
    double fH = model->H[0].dfunc(u, 0.0, /*log*/false);
    double age_hp = u + model->t0;
    double psi_u = calc_psi(age_hp, model);
    return fH * psi_u;
  }
  
  double indolent_inflow(double L_age, double U_age, Model* model) {
    
    if (!std::isfinite(U_age) || !std::isfinite(L_age) || U_age <= L_age)
      return 0.0;
    
    // convert to time since onset
    double L = L_age - model->t0;
    double U = U_age - model->t0;
    
    double width = U - L;
    
    auto integrand = [&](double x_unit) {
      // map from [0,1] -> [L,U]
      double u = L + width * x_unit;
      double jac = width;
      
      return indolent_kernel(u, model) * jac;
    };
    
    return gauss_legendre_32(integrand);
  }
  
  double integrand_prog_surv(double u, double last, Model* model) {
    
    double fH = model->H[0].dfunc(u, 0.0, /*log*/false);
    
    double age_hp = model->t0 + u;
    double psi_u = calc_psi(age_hp, model);
    double SP = model->P[0].pfunc(last - u, age_hp, /*lower.tail*/false, /*log*/false);
    return fH * (1.0 - psi_u) * SP;
  }
  
  double integral_prog_surv(double L_age, double U_age, double last, Model* model) {
    
    if (!std::isfinite(U_age) || !std::isfinite(L_age) || U_age <= L_age)
      return 0.0;
    
    // convert to time since onset
    double L = L_age - model->t0;
    double U = U_age - model->t0;
    
    double width = U - L;
    
    auto integrand = [&](double x_unit) {
      // map from [0,1] -> [L,U]
      double u = L + width * x_unit;
      double jac = width;
      
      return integrand_prog_surv(u, last - model->t0, model) * jac;
    };
    
    return gauss_legendre_32(integrand);
  }
  
  
  
  double D_j_A(Model* model, const std::vector<double>& time_points, int n_screens) {
    
    double A = 0.0;
    
    double beta = model->beta[0].parameterModel->parameters[0]->value;
    
    for (int k = 1; k <= n_screens; ++k) {
      
      double tmp = beta * std::pow(1.0 - beta, n_screens - k);
      
      double inflow = indolent_inflow(time_points[k-1], time_points[k], model);
      
      A += tmp * inflow;
    }
    
    return A;
  }
  
  double D_j_B(Model* model, std::vector<double>& time_points, int n_screens) {
    
    double B = 0.0;
    
    double beta = model->beta[0].parameterModel->parameters[0]->value;
    for (int k = 1; k <= n_screens; ++k) {
      double tmp = beta * std::pow(1.0 - beta, n_screens - k);
      B += tmp * integral_prog_surv(time_points[k-1], time_points[k], time_points[n_screens], model);
    }
    
    return B;
  }  
  
  /* circumvent computationally intensive integrals with parameters use
   * constant models */
  
  double prog_odx_fast(
      double age_min,
      double s,
      const std::vector<double>& lambda_ages,
      const std::vector<double>& lambda_rates,
      Model* model,
      double T_max)
  {
    const double S_screen = other_surv(age_min, s, lambda_ages, lambda_rates);
    
    auto integrand = [&](double y){
      const double t = y * T_max;
      return other_surv(age_min, s + t, lambda_ages, lambda_rates)
        * model->P[0].dfunc(t, 0.0, false)
        * T_max;
    };
    
    const double tail = gauss_legendre_32(integrand);
    
    return clamp01(S_screen - tail);
  }
  
  /*
   * Estimates the marginal probability that a progressive case detected
   * at a given screening age would die of other causes before becoming clinical.
   *
   * This function evaluates, for each scheduled screening age, the expected fraction
   * of progressive (non-indolent) cancers that would not present clinically
   * over an individual’s remaining lifetime due to death from other causes.
   *
   * The calculation integrates across the joint space of:
   *   - Age at transition from Healthy -> Preclinical (age_hp)
   *   - Preclinical sojourn duration (t)
   *   - Age-specific other-cause survival function
   *
   * The result represents:
   *    P(progressive, detected at screen s, dies before clinical onset),
   * averaged over the onset-age distribution and conditional on survival to the screen.
   *
   * Integration details:
   *  - Outer integral: over onset ages in [0, s - t0]
   *  - Inner integral: over preclinical durations t in [0, T_max]
   *  - Quadrature: 32-point Gauss–Legendre on both domains
   *
   * @param age_min Minimum age (typically first screening age).
   * @param screen_ages Vector of screening ages (in years).
   * @param lambda_ages Vector of ages corresponding to other-cause hazard rates.
   * @param lambda_rates Vector of other-cause hazard rates.
   * @param H_ddist Density function of the Healthy compartment.
   * @param H_params Parameter list for the Healthy compartment distribution.
   * @param psi_eval R function evaluating indolence probability psi(age_hp).
   * @param psi_params Parameters for the psi model.
   * @param P1_ddist Density of the primary preclinical (or fast/aggressive) distribution.
   * @param P1_pdist CDF of the primary preclinical distribution.
   * @param P1_params Parameters for the primary preclinical model.
   * @param t0 Risk onset age.
   * @param T_max Upper bound for integration over preclinical durations.
   * @param outer_power Power-transform tuning parameter for outer integration density.
   *
   * @return NumericVector of probabilities for each screening age.
   */
  NumericVector prog_odx_general(
      double age_min,
      NumericVector& screen_ages,
      std::vector<double>& lambda_ages,
      std::vector<double>& lambda_rates,
      Model* model,
      double T_max = 500.0,
      double outer_power = 1.0) {
    
    NumericVector out(screen_ages.size());
    
    // transformation factor controlling outer quadrature density
    const double beta = std::max(1e-6, outer_power);
    const double inv_beta = 1.0 / beta;
    
    // loop over all screening ages
    for (int i = 0; i < screen_ages.size(); ++i) {
      Rcpp::checkUserInterrupt();
      const double s = screen_ages[i]; // screening age
      const double U = s - model->t0; // time since risk onset
      if (!(U > 0)) { out[i] = 0.0; continue; }
      
      double Z = 0.0; // denominator (total kernel mass)
      double Num = 0.0; // numerator (weighted by death-before-clinical probability)
      
      // OUTER integral over onset age
      auto outer_integrand = [&](double y) -> double {
        // map y in [0,1] to u in [0,U] using power transform
        const double y_pow = std::pow(y, inv_beta);
        const double u = U * y_pow;
        const double jac_u = (U * inv_beta) * std::max(0.0, std::pow(y, inv_beta - 1.0));
        
        // age at disease onset
        const double age_hp = model->t0 + u;
        
        // f_H(u) density of healthy -> preclinical transition
        double fH = model->H[0].dfunc(u, 0.0, /*log*/false);
        if (!(fH > 0.0) || !R_finite(fH)) return 0.0;
        
        // (1 - psi) probability that transition is to progressive disease
        double psi_tau = calc_psi(age_hp, model);
        const double one_m_psi = 1.0 - psi_tau;
        if (!(one_m_psi > 0.0)) return 0.0;
        
        // remaining time to screen
        const double x = U - u;
        if (!(x > 0.0)) return 0.0;
        
        // survival in the preclinical state up to the screen
        double SP = model->P[0].pfunc(x, age_hp, /*lower.tail*/false, /*log*/false);
        
        // joint kernel: density of onset x non-indolence x preclinical survival
        const double kernel = fH * one_m_psi * clamp01(SP);
        if (!(kernel > 0.0) || !R_finite(kernel)) return 0.0;
        
        const double Soth_screen = other_surv(age_min, s, lambda_ages, lambda_rates);
        
        // INNER integral over sojourn time
        auto inner_integrand = [&](double y2) -> double {
          // map y2 in [0,1] to t in [0, T_max]
          const double t = y2 * T_max;
          const double jac_t = T_max;
          
          // other-cause survival to screen and after potential clinical onset
          const double Soth_after  = other_surv(age_min, s + t, lambda_ages, lambda_rates);
          // probability of death from other causes before clinical onset
          const double death_before_clinical = std::max(0.0, Soth_screen - Soth_after);
          
          // f_P(t): density of preclinical durations
          double fP = model->P[0].dfunc(t, age_hp, /*log*/false);
          return death_before_clinical * fP * jac_t;
        };
        
        // expected probability of death before clinical onset, given onset at u
        const double E_death_before_clinical = gauss_legendre_32(inner_integrand);
        // accumulate weighted totals
        Z += kernel * jac_u;
        Num += kernel * E_death_before_clinical * jac_u;
        return 0.0; // value not used; accumulating Z, Num
      };
      
      // evaluate outer integral
      (void)gauss_legendre_32(outer_integrand);
      
      // final probability = E[death before clinical | progressive onset]
      const double odx_prog = (Z > 1e-300) ? (Num / Z) : 0.0;
      out[i] = clamp01(odx_prog);
    }
    
    return out;
  }
  
  
  enum class TailPolicy {
    CAP,
    ZERO,
    EXTEND
  };
  
  inline double surv_fun_cpp(double age_min,
                             double t,
                             const std::vector<double>& ages,
                             const std::vector<double>& rates,
                             double terminal_age = std::numeric_limits<double>::infinity(),
                             TailPolicy tail = TailPolicy::EXTEND) {
    
    if (!std::isfinite(age_min) || !std::isfinite(t) || t <= age_min)
      return 1.0;
    
    // terminal-age policy
    if (std::isfinite(terminal_age)) {
      if (tail == TailPolicy::ZERO && t > terminal_age)
        return 0.0;
      if (tail == TailPolicy::CAP && t > terminal_age)
        t = terminal_age;
      // EXTEND does nothing
    }
    
    const double max_age = ages.back();
    const double t_tab = std::min(t, max_age);
    
    double base_hh = 0.0;
    double left = age_min;
    
    while (left < t_tab) {
      
      // find interval index j such that ages[j] <= left < ages[j+1]
      auto it = std::upper_bound(ages.begin(), ages.end(), left);
      int idx = std::max(0, int(it - ages.begin()) - 1);
      
      double right = t_tab;
      if (idx + 1 < (int)ages.size())
        right = std::min(right, ages[idx + 1]);
      
      if (right > left)
        base_hh += rates[idx] * (right - left);
      
      left = right;
    }
    
    // tail beyond last tabulated age (extend last hazard)
    if (t > max_age) {
      double last_rate = rates.back();
      base_hh += last_rate * (t - max_age);
    }
    
    return std::exp(-base_hh);
  }
}

// [[Rcpp::export]]
List predict_ODX(List H, List P, List psi, List beta, double t0, 
                 Rcpp::NumericVector screening_schedule_r,
                 Rcpp::NumericVector lambda_ages_r,
                 Rcpp::NumericVector lambda_rates_r,
                 Rcpp::NumericMatrix posteriors,
                 Rcpp::CharacterVector posterior_column_names) {
  
  Model model = Model();
  try {
    model = Model(H, P, psi, beta, t0);
  } catch (const NegativeRateError& error) {
    stop("could not create model");
  }
  
  const int J = screening_schedule_r.size();
  const int M = posteriors.nrow();
  if (posteriors.ncol() != posterior_column_names.size())
    stop("posterior_column_names length must match posteriors.ncol()");
  
  std::vector<double> screens(J);
  for (int j = 0; j < J; ++j) screens[j] = screening_schedule_r[j];
  
  std::vector<int> ord(lambda_ages_r.size());
  for (int i = 0; i < (int)ord.size(); ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b){ return lambda_ages_r[a] < lambda_ages_r[b]; });
  
  std::vector<double> haz_ages(ord.size()), haz_rates(ord.size());
  for (int i = 0; i < (int)ord.size(); ++i) {
    haz_ages[i]  = lambda_ages_r[ord[i]];
    double r     = lambda_rates_r[ord[i]];
    haz_rates[i] = (std::isfinite(r) && r > 0.0) ? r : 0.0;
  }
  
  std::unordered_map<std::string,int> col_index;
  col_index.reserve(posterior_column_names.size() * 2);
  
  for (int c = 0; c < posterior_column_names.size(); ++c) {
    col_index.emplace(Rcpp::as<std::string>(posterior_column_names[c]), c);
  }
  
  const int Pn = (int)model.parameters.size();
  std::vector<int> pcol(Pn, -1);
  
  for (int i = 0; i < Pn; ++i) {
    const std::string tag = model.parameters[i]->tag;
    auto it = col_index.find(tag);
    if (it == col_index.end()) {
      stop(("posterior matrix missing required parameter column: " + tag).c_str());
    }
    pcol[i] = it->second;
  }
  
  NumericMatrix screen_total(J, M);
  NumericMatrix screen_ind(J, M);
  NumericMatrix screen_prog(J, M);
  
  NumericMatrix SD_stat(J, M);
  NumericMatrix ODX_ind_stat(J, M);
  NumericMatrix ODX_prog_stat(J, M);
  
  NumericVector hh1(J);
  for (int j = 0; j < J; ++j) {
    hh1[j] = surv_fun_cpp(screens[0], screens[j], haz_ages, haz_rates,
                          /*terminal_age=*/std::numeric_limits<double>::infinity(),
                          TailPolicy::EXTEND);
  }
  
  std::vector<double> time_points;
  time_points.reserve(J + 1);
  
  for (int m = 0; m < M; ++m) {
    Rcpp::checkUserInterrupt();
    for (int k = 0; k < Pn; ++k) {
      model.parameters[k]->value = posteriors(m, pcol[k]);
    }

    NumericVector AFS(1, screens[0]);
    const NumericVector log_norm =
      normalization::compute_cp_log_generic(AFS, model.t0, model.H[0], model.P[0], model.psi[0]);
    const double norm_inv = std::exp(-log_norm[0]);

    time_points.clear();
    time_points.push_back(model.t0);
    
    for (int ns = 0; ns < J; ++ns) {
      time_points.push_back(screens[ns]);
      
      const int n_screens = ns + 1;
      
      const double dja = D_j_A(&model, time_points, n_screens);
      const double djb = D_j_B(&model, time_points, n_screens);
      
      const double tot = (dja + djb) * norm_inv;
      
      screen_total(ns, m) = tot;
      screen_ind(ns, m)   = dja * norm_inv;
      screen_prog(ns, m)  = djb * norm_inv;
    }
    
    NumericVector prog_vec(J);
    if (model.P[0].parameterModel->typeid_ == Pkg::Model::CONSTANT &&
        model.psi[0].parameterModel->typeid_ == Pkg::Model::CONSTANT) {
      for (int ns = 0; ns < J; ++ns) {
        prog_vec[ns] = prog_odx_fast(screens[0], screens[ns], 
                                     haz_ages, haz_rates, &model, 
                                     /*T_max*/500);
      }
    } else {
      prog_vec = prog_odx_general(
        screens[0], screening_schedule_r,
        /*lambda_ages*/ haz_ages,
        /*lambda_rates*/ haz_rates,
        &model,
        /*T_max*/ 500.0,
        /*outer_power*/ 1.0);
    }
    
      
    for (int j = 0; j < J; ++j) {
      double v = prog_vec[j];
      if (!R_finite(v) || v < 0.0) v = 0.0;
      prog_vec[j] = v;
    }
      
    for (int j = 0; j < J; ++j) {
      SD_stat(j, m) = screen_total(j, m) * hh1[j];
      ODX_ind_stat(j, m) = screen_ind(j, m) * hh1[j];
      ODX_prog_stat(j, m) = screen_prog(j, m) * prog_vec[j];
    }
  }
  
  return List::create(
    Named("ODX_ind_stat")  = ODX_ind_stat,
    Named("ODX_prog_stat") = ODX_prog_stat,
    Named("SD_stat")       = SD_stat
  );
}