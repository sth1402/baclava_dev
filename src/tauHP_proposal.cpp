#include <Rcpp.h>
#include "Model.h"
#include "tauHP_proposal.h"
#include "NegativeRateError.h"
#include "EPGroup_State.h"
#include "StateParameter.h"
#include "Pkg_types.h"

using namespace Rcpp;

// create the weighting matrix once
NumericVector compute_prob_tau_prep(const EPGroup_State& gs, const Model* model) {
  
  NumericVector wgt(gs.n_intvl, 1.0);
  
  // Screening adjustment (if applicable)
  if (gs.age_screen.size() == 0) return wgt;
  
  // assign beta to each assessment
  
  NumericVector beta_by_type(model->beta.size());
  for (int i = 0; i < model->beta.size(); ++i) {
    beta_by_type[i] = model->beta[i].parameterModel->par_from_tau(0.0);
  }
  
  if (gs.endpoint_type == Pkg::EP::PRECLINICAL) {
    // weighting is (1-b)^{m-1}b, ..., (1-b)b, b
    int start = 0;
    for (int i = 0; i < gs.n; ++i) {
      
      // if only 1 interval set to 1
      if (gs.ep_values[i].size() == 2) {
        start++;
        continue;
      }
      
      int end = start + gs.ep_values[i].size() - 1;
      double res;
      
      int ns = gs.screen_types[i].size() - 1;
      res = beta_by_type[gs.screen_types[i][ns]];
      ns--;
      for (int j = end - 1; j > start; --j) {
        wgt[j] *= res;
        res *= (1.0 - beta_by_type[gs.screen_types[i][ns]]);
        ns--;
      }
      wgt[start] *= res;
      start = end;
    }
  } else if (gs.endpoint_type == Pkg::EP::CLINICAL) {
    // weighting is (1-b)^{m}, ..., (1-b), 1
    int start = 0;
    for (int i = 0; i < gs.n; ++i) {
      // last interval for clinical cases is (last_screen -> diagnosis)
      // this interval is not weighted by beta
      
      int end = start + gs.ep_values[i].size() - 1;
      int ns = gs.screen_types[i].size() - 1;
      double res = (1.0 - beta_by_type[gs.screen_types[i][ns]]);
      ns--;
      for (int j = end - 2; j > start; --j) {
        wgt[j] *= res;
        res *= (1.0 - beta_by_type[gs.screen_types[i][ns]]);
        ns--;
      }
      wgt[start] *= res;
      start = end;
    }
  } else if (gs.endpoint_type == Pkg::EP::CENSORED) {
    // 3 versions
    // 1: weighting is (1-b)^{m+1, }..., (1-b), 1
    // 2: weighting is 1, (1-b)^{m-1}, ..., (1-b), 1
    // 3: weighting is (1-b)^{m}, ..., (1-b), 1, 1
    int start = 0;

    for (int i = 0; i < gs.n; ++i) {
      // last two intervals for censored cases are (last_screen -> exit) (exit -> Inf)
      // these are not weighted by beta
      
      int end = start + gs.ep_values[i].size() - 1;
      int ns = gs.screen_types[i].size() - 1;
      double res = (1.0 - beta_by_type[gs.screen_types[i][ns]]);
      
      ns--;
      for (int j = end - 3; j > start; --j) {
        wgt[j] *= res;
        res *= (1.0 - beta_by_type[gs.screen_types[i][ns]]);
        ns--;
      }
      wgt[start] *= res;
      start = end;
    }
  }
  
  return wgt;
  
}

NumericVector compute_prob_tau(
    const EPGroup_State& gs, 
    const NumericVector& tau, 
    const NumericVector& wgt, 
    const Model* model) {

  // Evaluate interval-wise probabilities
  
  NumericVector prob_tau;
  if (gs.endpoint_type != Pkg::EP::CLINICAL) {
    prob_tau = model->H[0].pfunc_ab_optimized(gs, NumericVector(gs.n, 0.0), model->t0, wgt);
  } else {
    prob_tau = model->P[gs.iP].pfunc_ab_optimized(gs, tau, model->t0, wgt);
  }

  if (is_true(any(prob_tau < 0))) 
      throw NegativeRateError("probability is negative");

  return prob_tau;

}

static inline int draw_discrete_01(const double* p, int K, double u) {
  double cum = 0.0;
  for (int k = 0; k < K; ++k) {
    cum += p[k];
    if (u <= cum) return k;
  }
  return K - 1;
}


/*
 * Samples new proposed values of tau (the latent age at which
 *   a detectable tumor appears) for each subject, based on the interval-specific 
 *   proposal probabilities constructed in compute_prob_tau_obj().
 *   
 * This function forms the proposal kernel for the tau Metropolis-Hastings 
 *   update, and incorporates subject-specific screening data and endpoint type.
 *   
 * gs: A EPGroup_State of subject-specific data, including endpoint time, 
 *   possible intervals for tau, screening history, etc.
 * tau: Vector of current tau transition ages.
 * prob_tau: The discrete distribution over tau intervals for each subject 
 *   (output of compute_prob_tau_obj()).
 * model: A pointer to the current model object, containing access to H/P 
 *   distributions and screen sensitivity models.
 */
NumericVector rprop_tau(const EPGroup_State& gs, 
                        const NumericVector& tau, 
                        const NumericVector& prob_tau,
                        const Model* model) {
  
  /*
   * Select interval for each subject
   *   Uses the subject-specific discrete probability vector from prob_tau to 
   *     draw an index (k_new[i]) indicating the chosen interval [a_i, b_i].
   *   Uses ep (age-at-time-of boundaries) to retrieve the corresponding lower 
   *     and upper bounds for each selected interval.
   */

  IntegerVector k_new(gs.n, -1);
  NumericVector vec_a = no_init(gs.n);
  NumericVector vec_b = no_init(gs.n);
  
  int start = 0;
  for (int i = 0; i < gs.n; ++i) {
    int end = start + gs.ep_values[i].size() - 1;
    int K = gs.ep_values[i].size() - 1;
    double u = R::runif(0.0, 1.0);
    k_new[i]  = draw_discrete_01(&(prob_tau[start]), K, u);
    
    vec_a[i] = gs.ep_values[i][k_new[i]];
    vec_b[i] = gs.ep_values[i][k_new[i] + 1];
    start = end;
  }

  if (is_true(any(vec_b < vec_a))) 
    stop("upper boundary below lower boundary; contact developer");
  
  NumericVector tau_new = no_init(gs.n);
  
  if (gs.endpoint_type == Pkg::EP::PRECLINICAL) {
    // sample tau from (interval_lower, interval_upper]
    NumericVector no_param(gs.n, 0.0);
    NumericVector sojourn_H_new = model->H[0].rfunc_trunc(vec_a - model->t0, 
                                                          vec_b - model->t0, 
                                                          no_param);
      
    tau_new = model->t0 + sojourn_H_new;
    
    if (is_true(any(tau_new >= gs.endpoint_times))) {
      tau_new = pmin(tau_new, gs.endpoint_times - 1e-12);
    }
    
  } else if (gs.endpoint_type == Pkg::EP::CENSORED) {
      
      if(gs.age_screen.size() != 0) {
        // Censor cases
        // sample tau from (interval_lower, interval_upper]
        NumericVector no_param(gs.n, 0.0);
        NumericVector sojourn_H_new = model->H[0].rfunc_trunc(vec_a - model->t0, 
                                                              vec_b - model->t0, 
                                                              no_param);
        
        tau_new = model->t0 + sojourn_H_new;
        // if tau sampled from [EP, Inf) set to infinite point mass
        tau_new[tau_new > gs.endpoint_times] = R_PosInf;
        
      } else {
      /* Censored cases that have no screens for these individuals, 
       * endpoint = (t0, T, Inf) for individuals with k_new = 1, the sampling 
       * interval is in the |T, inf) interval. These are unchanged. For 
       * individuals with k_new = 0, the sampling interval is in the
       * |t0, T) interval. These will now sample from an exponential
       */
      tau_new = rep(R_PosInf, gs.n);
      
      LogicalVector finite_intv = k_new == 0;
      if (is_true(any(finite_intv))) {

        NumericVector scale = model->P[gs.iP].compute_optimal_scale(tau[finite_intv]);
        NumericVector scale_clamped = pmax(scale, 1e-12);
        NumericVector vv(scale_clamped.size(), 0.0);
        std::transform(scale_clamped.begin(), 
                       scale_clamped.end(), 
                       vv.begin(), [=](double p){ return rexp(1, 1.0 / p)[0]; });
        
        NumericVector vec_T = gs.endpoint_times[finite_intv];
        NumericVector tau_new = vec_T - vv;
        
        // enforce tau in [t0, T)
        tau_new = pmax(model->t0, pmin(tau_new, vec_T - 1e-12));
        
        tau_new[finite_intv] = tau_new;
      }
    }

  } else if (gs.endpoint_type == Pkg::EP::CLINICAL) {
    // Clinical cases
    //sample tau_p from (interval_lower, interval_upper]
    NumericVector sojourn_P_new = model->P[gs.iP].rfunc_trunc(gs.endpoint_times - vec_b,
                                                              gs.endpoint_times - vec_a,
                                                              tau);
    // Guard against numerical underflow in truncated Weibull sampler
    // when gap between last screen and clinical diagnosis is very short (<0.01 years)
    // Biologically implausible sojourns of this length are floored at 1e-3
    sojourn_P_new = pmax(sojourn_P_new, NumericVector(sojourn_P_new.size(), 1e-3));
    
    tau_new = pmax(gs.endpoint_times - sojourn_P_new, model->t0 + 0.001);

  }

  if (is_true(any(!(is_finite(tau_new) | is_infinite(tau_new)))) ) {
      throw NegativeRateError("proposed tau is NA or NaN");
  }
  
  return tau_new;
}

/*
 * Computes the log-probability density under the proposal kernel of a given 
 *   tau vector, used to evaluate the Metropolis-Hastings proposal 
 *   ratio during the update of tau.
 *   
 * This function reflects the forward or reverse probability of transitioning 
 *   from one tau value to another under a proposal kernel based on:
 *     Discrete probabilities for interval selection
 *     Continuous densities within selected intervals (truncated Weibull)
 *     
 *  gs: EPGroup_State Individual-level data including endpoint time, interval 
 *    boundaries, screen info, etc.
 *  tau_eval: tau values whose probabilities are being evaluating under a proposal kernel
 *  prob_tau_eval: interval probs from kernel centered at new tau_eval
 *  tau_ref: reference (conditioning) transition time used by the proposal kernel.
 *  model: Current model pointer containing compartmental distributions and parameters.
 */
NumericVector dlog_prop_tau(const EPGroup_State& gs,
                            const NumericVector& tau_eval,
                            const NumericVector& prob_tau_eval,
                            const NumericVector& tau_ref,
                            const Model* model) {
  
  /*
   * Locate the proposal interval for each subject
   *   Identifies which interval [a,b] each tau[i] falls in.
   *   Uses this to:
   *     Extract the log-probability of selecting the interval (from prob_tau)
   *     Retrieve the interval boundaries [a,b]
   */
  
  NumericVector dlog_tau = no_init(gs.n);
  NumericVector dlog_k = no_init(gs.n);
  
  if (gs.age_screen.size() > 0 || 
      (gs.age_screen.size() == 0 && gs.endpoint_type == Pkg::EP::CLINICAL)) {
    
    // k_new is the last screen before the current estimated transition time
    IntegerVector k_new(gs.n, -1);
    NumericVector prob_vec = no_init(gs.n);
    NumericVector vec_a = no_init(gs.n);
    NumericVector vec_b = no_init(gs.n);

    int start = 0;
    for (int i = 0; i < gs.n; i++) {
      const auto& ep = gs.ep_values[i];
      
      int end = start + ep.size() - 1;
      
      if (std::isfinite(tau_eval[i])) {
        
        auto it = std::upper_bound(ep.begin(), ep.end(), tau_eval[i]);
        
        k_new[i] = int(it - ep.begin()) - 1;
        if (k_new[i] < 0) throw NegativeRateError("Invalid k_new");
      } else {
        k_new[i] = ep.size() - 2;
      }
      
      if (((k_new[i] + start) >= prob_tau_eval.size()) ||
          k_new[i] >= ep.size() || (k_new[i]+1) >= ep.size()) {
        Rcpp::Rcerr << "Out of bounds for individual " << i << "\n"
                    << "  endpoint_type: " << static_cast<int>(gs.endpoint_type) << "\n"
                    << "  k_new[i]: " << k_new[i] << "\n"
                    << "  ep.size(): " << ep.size() << "\n"
                    << "  prob_tau_eval.size(): " << prob_tau_eval.size() << "\n"
                    << "  start: " << start << "\n"
                    << "  k_new[i]+start: " << k_new[i]+start << "\n"
                    << "  tau_eval[i]: " << tau_eval[i] << "\n"
                    << "  ep[0]: " << ep[0] << "\n"
                    << "  ep[ep.size()-1]: " << ep[ep.size()-1] << "\n";
        stop("vector went out of bounds; contact developer");
      }
      prob_vec[i] = prob_tau_eval[k_new[i] + start];
      vec_a[i] = ep[k_new[i]] - model->t0;
      vec_b[i] = ep[k_new[i] + 1] - model->t0;
      start = end;
    }

     dlog_k = log(prob_vec);
    
    if (gs.endpoint_type == Pkg::EP::PRECLINICAL) {
      NumericVector no_param(gs.n, 0.0);
      dlog_tau = model->H[0].dfunc_trunc(tau_eval - model->t0, 
                                         vec_a, 
                                         vec_b, 
                                         no_param,
                                         /*uselog=*/true);
    } else if (gs.endpoint_type == Pkg::EP::CENSORED) {
      NumericVector no_param(gs.n, 0.0);
      dlog_tau = model->H[0].dfunc_trunc(tau_eval - model->t0, 
                                         vec_a, 
                                         vec_b, 
                                         no_param,
                                         /*uselog=*/true);
      dlog_tau[is_infinite(tau_eval)] = 0.0;
    } else if (gs.endpoint_type == Pkg::EP::CLINICAL) {

      NumericVector Tp = gs.endpoint_times - tau_eval;
      NumericVector a = gs.endpoint_times - model->t0 - vec_b;
      NumericVector b = gs.endpoint_times - model->t0 - vec_a;

      dlog_tau = model->P[gs.iP].dfunc_trunc(Tp, a, b, tau_ref, /*uselog=*/true);

    }
  } else if (gs.age_screen.size() == 0  && gs.endpoint_type != Pkg::EP::CLINICAL) {
    /*
     * handles:
     * A degenerate two-interval scenario: [t0, T] and [T, Inf)
     *  Samples from exponential within [t0, T]
     *  Assumes point-mass at Inf otherwise (density = 0; log-density = 0)
     */
    NumericVector tmp_vec(gs.n);
    int start = 0;
    for (int i = 0; i < gs.n; ++i) {
      int end = start + gs.ep_values[i].size() - 1;
      dlog_tau[i] = 0.0;
      tmp_vec[i] = prob_tau_eval[end-1];
      start = end;
    }
    
    LogicalVector tst = tau_eval < gs.endpoint_times;
    if (is_true(any(tst))) {
      NumericVector scale = model->P[gs.iP].compute_optimal_scale(tau_ref[tst]);
      
      start = 0;
      for (int i = 0; i < gs.n; ++i) {
        int end = start + gs.ep_values[i].size() - 1;
        if (tst[i]) {
          double vec_a = gs.endpoint_times[i] - tau_eval[i];
          dlog_tau[i] = R::dexp(vec_a, 1.0 / scale[i], true);
          tmp_vec[i] = prob_tau_eval[start];
        }
        start = end;
      }
    }    
    dlog_k = log(tmp_vec);
  } else {
    stop("unrecognized endpoint type");
  }
  
  return dlog_k + dlog_tau;
}