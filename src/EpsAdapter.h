#ifndef EPSADAPTER_H
#define EPSADAPTER_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "Verbosity.h"

/*
 * EpsAdapter
 * ----------
 *
 * Adaptive step-size (global scale) controller for random-walk Metropolis
 * and HMC-like samplers.
 *
 * This class adapts a *single global step size* `s = exp(log_eps)` that
 * scales proposal jumps independently of covariance geometry. It is intended
 * to be used during warmup only and frozen before sampling.
 *
 * The adapter operates entirely on the log scale to ensure positivity
 * and numerical stability.
 *
 * ---------------------------------------------------------------------------
 * Conceptual role
 * ---------------------------------------------------------------------------
 *
 * In a block-wise proposal of the form
 *
 *     x_new = x + s * (D * L) * z ,   z ~ N(0, I)
 *
 * this adapter controls `s`, while `CovAdapter` controls the geometry
 * (D, L). These responsibilities are deliberately separated:
 *
 *   - CovAdapter learns *shape*
 *   - EpsAdapter learns *overall scale*
 *
 * ---------------------------------------------------------------------------
 * Adaptation modes
 * ---------------------------------------------------------------------------
 *
 * Two adaptation schemes are supported:
 *
 * 1. "fast" mode (default)
 *    ---------------------
 *    Dual-averaging (Nesterov-style) adaptation as used in Stan
 *    (Hoffman & Gelman, 2014).
 *
 *    This method:
 *      - Targets a desired acceptance probability `delta`
 *      - Uses an accumulated error term (H_accum)
 *      - Adapts aggressively early in warmup
 *      - Smoothly stabilizes via a decaying learning rate
 *
 *    The update equations are:
 *
 *      H_t      = (1 - rho_t) * H_{t-1} + rho_t * (delta - alpha_bar)
 *      log s_t  = mu - (sqrt(t) / gamma) * H_t
 *      log sbar_t = t^{-kappa} * log s_t + (1 - t^{-kappa}) * log sbar_{t-1}
 *
 *    where:
 *      - alpha_bar is the batch-averaged acceptance probability
 *      - mu is an anchoring location (default: log(10 * s0))
 *      - gamma controls aggressiveness
 *      - kappa controls smoothing
 *      - t0 delays early adaptation for stability
 *
 * 2. "slow" mode
 *    ------------
 *    Robbins–Monro stochastic approximation on log step size.
 *
 *    This method:
 *      - Updates more conservatively
 *      - Uses diminishing step sizes
 *      - Is appropriate when acceptance noise is high
 *
 *    The update is:
 *
 *      log s_t = log s_{t-1} + rho_t (alpha_bar − delta)
 *
 *    with exponential smoothing applied to form log_eps_bar.
 *
 * ---------------------------------------------------------------------------
 * Acceptance handling
 * ---------------------------------------------------------------------------
 *
 * - Adaptation is based on *observed acceptance probability*
 * - Acceptance values are clipped to [0, 1]
 * - Updates are performed in small batches to reduce noise
 * - Adaptation does NOT depend on proposal acceptance/rejection alone,
 *   but on the realized acceptance probability
 *
 * ---------------------------------------------------------------------------
 * Freezing behavior
 * ---------------------------------------------------------------------------
 *
 * When `freeze()` is called:
 *
 *   - Adaptation stops permanently
 *   - The final step size is taken as `exp(log_eps_bar)`
 *   - Bounds [s_floor, s_ceil] are enforced
 *
 * ---------------------------------------------------------------------------
 * Defaults and dimensional heuristics
 * ---------------------------------------------------------------------------
 *
 * If `delta` is not specified, it is chosen based on block dimension:
 *
 *   d = 1     delta = 0.44
 *   d <= 4    delta = 0.35
 *   d > 4     delta = 0.234
 *
 * ---------------------------------------------------------------------------
 * Design principles
 * ---------------------------------------------------------------------------
 *
 * - Log-scale adaptation for numerical safety
 * - Geometry-independent global scaling
 * - Batch-averaged acceptance to reduce jitter
 * - Smooth finalization to avoid warmup artifacts
 */

struct EpsAdapter {
  
  EpsAdapter() {};
  
  /*
   * Create from provided R based list
   *   
   *   initial_eps: s0 from block definition
   *   adaptive list: delta, warmup, t0, kappa, s_floor, s_ceil
   *   p: number of parameters in update block
   */
  EpsAdapter(const double initial_eps, 
             const Rcpp::List& adaptive_list, 
             const int p) : EpsAdapter() {
    
    Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),
                  "[epsadapt] initializing from list");
    
    if (adaptive_list.size() == 0) {
      this->frozen = true;
      this->log_eps = std::min(std::max(std::log(initial_eps), 
                                        std::log(this->s_floor)), 
                                        std::log(this->s_ceil));
      this->log_eps_bar = log_eps;
      return;
    }

    if (!adaptive_list.containsElementNamed("mode")) {
      this->mode = "fast";
    } else {
      this->mode = Rcpp::as<std::string>(adaptive_list["mode"]);
    }
    if (this->mode != "fast" && this->mode != "slow") Rcpp::stop("mode must be 'fast' or 'slow'");
    
    // adaptation controls
    if (!adaptive_list.containsElementNamed("delta")) {
      // if delta not specified, use defaults based on p
      this->delta = set_dim_target_(p);
    } else {
      // if delta specified by user, keep it.
      this->delta = Rcpp::as<double>(adaptive_list["delta"]); 
      if (Rcpp::NumericVector::is_na(this->delta) || std::isnan(this->delta) )
        this->delta = set_dim_target_(p);
    }

    if (this->delta <= 0.0 || this->delta >= 1.0) Rcpp::stop("delta must be in (0,1)");
    
    if (adaptive_list.containsElementNamed("kappa")) 
      this->kappa = Rcpp::as<double>(adaptive_list["kappa"]);
    if (this->kappa > 1.0 || this->kappa <= 0.5) Rcpp::stop("kappa must be in (0.5, 1.0]");
    if (adaptive_list.containsElementNamed("t0")) 
      this->t0 = Rcpp::as<double>(adaptive_list["t0"]);
    if (this->t0 <= 0.0) Rcpp::stop("t0 must be positive");

    // terms for fast adaptation    
    if (adaptive_list.containsElementNamed("gamma")) 
      this->gamma = Rcpp::as<double>(adaptive_list["gamma"]);
    if (this->mode == "fast" && this->gamma <= 0.0) Rcpp::stop("gamma must be positive");
    if (adaptive_list.containsElementNamed("mu")) 
      this->mu = Rcpp::as<double>(adaptive_list["mu"]);
    
    // adaptation procedural controls
    
    if (adaptive_list.containsElementNamed("s_floor")) 
      this->s_floor = Rcpp::as<double>(adaptive_list["s_floor"]);
    if (this->s_floor <= 0.0) Rcpp::stop("s_floor cannot be negative");
    if (adaptive_list.containsElementNamed("s_ceil")) 
      this->s_ceil  = Rcpp::as<double>(adaptive_list["s_ceil"]);
    if (this->s_ceil <= this->s_floor) Rcpp::stop("s_ceil cannot be <= s_floor");
    if (adaptive_list.containsElementNamed("batch")) 
      this->batch  = Rcpp::as<int>(adaptive_list["batch"]);
    if (this->batch < 1) Rcpp::stop("batch must be >= 1");
    
    if (initial_eps <= 0.0) Rcpp::stop("s0 must be positive");
    this->log_eps = std::min(std::max(std::log(initial_eps), 
                                      std::log(this->s_floor)), 
                                      std::log(this->s_ceil));
    this->log_eps_bar = log_eps;
    this->it = 0;
    if (R_IsNA(this->mu) || !std::isfinite(this->mu)) this->mu = std::log(10.0 * initial_eps);
    
    std::string msg = "[epsadapt] delta=" + std::to_string(this->delta)
      + " kappa=" + std::to_string(this->kappa)
      + " t0=" + std::to_string(this->t0);
      
      if (this->mode == "fast") {
        msg = msg + " mu=" + std::to_string(this->mu)
        + " gamma=" + std::to_string(this->gamma);
        
      }
      msg = msg + 
        " s_bound=(" + std::to_string(this->s_floor) + "," + std::to_string(this->s_ceil) + ")"
      + " batch=" + std::to_string(this->batch);
      
      Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),  msg);
      
  }
  
  EpsAdapter(const double initial_eps, const int p) {
    
    Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),
                  "[epsadapt] initializing without list");
    
    this->delta = set_dim_target_(p);

    std::string msg = "[epsadapt] delta=" + std::to_string(this->delta)
      + " kappa=" + std::to_string(this->kappa)
      + " t0=" + std::to_string(this->t0)
      + " s_bound=(" + std::to_string(this->s_floor) + "," + std::to_string(this->s_ceil) + ")"
      + " batch=" + std::to_string(this->batch);
      
      Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),  msg);
      
    this->log_eps = std::log(initial_eps);
    this->log_eps_bar = log_eps;
    this->it = 0;
  }
  
  double update(double alpha_hat) {
    if (this->frozen) return(std::exp(this->log_eps));
    if (this->mode == "fast") {
      return update_fast(alpha_hat);
    } else {
      return update_slow(alpha_hat);
    }
  }
  
  double current() const { return std::exp(this->log_eps); }
  double current_bar() const { return std::exp(this->log_eps_bar); };
  double getDelta() const { return this->delta; }
  double freeze() { 
    this->frozen = true;
    return this->finalize(); 
  }
  bool isFrozen() { return this->frozen; }

private:
  
  bool frozen = false;
  
  std::string mode;
  
  double delta = 0.44; // target acceptance
  double t0 = 10.0; // stability; larger values essentially delay adaptation
  double kappa = 0.75; // decay exponent
  
  // for more aggressive updating
  double H_accum = 0.0; // integral feedback term
  double gamma = 0.2; // tuneable gain/aggressiveness larger -> faster/noisier
  double mu = R_NaReal; // anchor
  
  // bound s
  double s_floor = 1e-4;
  double s_ceil = 2.5;
  
  int batch = 1;
  int acc_cnt = 0; // batch counter
  double acc_sum = 0.0; // batch accumulator
  
  // state
  int it = 0;
  double log_eps = 0.0; // current log step
  double log_eps_bar = 0.0; // smoothed log step
  
  // delta default depends on number of parameters in update block
  static double set_dim_target_(int d) {
    if (d == 1) return 0.44;
    if (d <= 4) return 0.35;
    return 0.234;
  }
  
  // Robbins–Monro stochastic approximation for log-step-size adaptation
  double update_slow(double alpha_hat) {
    // accumulate acceptance in small batches to reduce jitter
    double a = alpha_hat;
    if (a < 0.0) a = 0.0;
    if (a > 1.0) a = 1.0;
    
    acc_sum += a;
    acc_cnt += 1;
    
    // only adapt when a full batch is collected
    if (acc_cnt >= batch) {
      ++it;
      
      const double a_bar = acc_sum / static_cast<double>(acc_cnt);
      acc_sum = 0.0;
      acc_cnt = 0;
      
      // Robbins–Monro
      double rho_t = 1.0 / (static_cast<double>(this->it) + this->t0);
      this->log_eps += rho_t * (a_bar - this->delta);
      
      // clamp in log-space
      const double log_floor = std::log(this->s_floor);
      const double log_ceil  = std::log(this->s_ceil);
      this->log_eps = std::min(std::max(this->log_eps, log_floor), log_ceil);
      
      double t_k = std::pow(static_cast<double>(this->it), -this->kappa);
      this->log_eps_bar = t_k * this->log_eps + (1.0 - t_k) * this->log_eps_bar;
      
      Verbosity::instance().log(TRACE, "[epsadapt] s=", std::exp(this->log_eps));
      
    } else {
      Verbosity::instance().log(TRACE, "[epsadapt] s unchanged");
    }
    
    return std::exp(this->log_eps);
  }
  
  // Dual averaging (Nesterov) for step-size adaptation (Hoffman & Gelman; Stan).
  double update_fast(double alpha_hat) {
    // clip
    double a = alpha_hat;
    if (a < 0.0) a = 0.0;
    if (a > 1.0) a = 1.0;
    
    // accumulate within batch
    acc_sum += a;
    acc_cnt += 1;
    
    if (acc_cnt >= batch) {
      ++it;
      const double a_bar = acc_sum / static_cast<double>(acc_cnt);
      acc_sum = 0.0;
      acc_cnt = 0;
      
      double rho_t = 1.0 / (static_cast<double>(this->it) + this->t0);
      this->H_accum = (1.0 - rho_t) * this->H_accum + (this->delta - a_bar) * rho_t;
      
      this->log_eps = this->mu - (std::sqrt(static_cast<double>(this->it)) / this->gamma) * this->H_accum;
      
      const double log_floor = std::log(this->s_floor);
      const double log_ceil  = std::log(this->s_ceil);
      this->log_eps = std::min(std::max(this->log_eps, log_floor), log_ceil);
      
      double t_k = std::pow(static_cast<double>(this->it), -this->kappa);
      this->log_eps_bar = t_k * this->log_eps + (1.0 - t_k) * this->log_eps_bar;
      
      Verbosity::instance().log(TRACE, "[epsadapt-fast] s=", std::exp(this->log_eps));
    } else {
      Verbosity::instance().log(TRACE, "[epsadapt-fast] s unchanged");
    }
    
    return std::exp(this->log_eps);
  }

  double finalize() const {
    double frozen_log = this->log_eps_bar;
    const double log_floor = std::log(s_floor);
    const double log_ceil  = std::log(s_ceil);
    if (frozen_log < log_floor) frozen_log = log_floor;
    if (frozen_log > log_ceil)  frozen_log = log_ceil;
    
    Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),  
                        "[epsadapt] final s=", std::exp(frozen_log));
    
    return std::exp(frozen_log);
  }
  
};
#endif
