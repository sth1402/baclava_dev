#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include "CorrBlock.h"
#include "CovAdapter.h"
#include "EpsAdapter.h"
#include "Proposal_State.h"
#include "StateParameter.h"
#include "Verbosity.h"

using namespace Rcpp;

// tracking data for debugging
void CorrBlock::adapter_track(double a_hat) {
  const double s_now = eps_adapter.current();
  const double s_bar = eps_adapter.current_bar();
  s_evolution.push_back(s_now);
  sbar_evolution.push_back(s_bar);
  
  ahat_evolution.push_back(a_hat);
  
  if (this->learn_D) 
    d_evolution.emplace_back(Rcpp::as<std::vector<double>>(this->D));
  
  if (this->learn_L) {
    // L diagnostics
    // upper triangle of L (off-diagonal evolution)
    
    int n = this->L.nrow();
    int k = n * (n - 1) / 2;   // number of elements above diagonal
    std::vector<double> off_diags;
    off_diags.reserve(k);
    
    for (int j = 1; j < n; ++j) {
      for (int i = 0; i < j; ++i) {
        off_diags.push_back(this->L(i, j));
      }
    }
    l_evolution.emplace_back(off_diags);
  }
}

// proposal
NumericVector CorrBlock::propose(const NumericVector& w_old_block,
                                 const Proposal_State& params) const {
  const int d = static_cast<int>(w_old_block.size());

  NumericVector w_new(d);
  if (d == 0) return w_new;
  
  // z ~ N(0, I_d)
  // runif option is included only to match the proposal step of the
  // original implementation - it should not be used in joint update steps
  NumericVector z(d);
  for (int j = 0; j < d; ++j) {
    if (params.use_rnorm[j]) {
      z[j] = R::rnorm(0.0, 1.0);
    } else if (params.use_runif[j]) {
      z[j] = R::runif(-1.0, 1.0);
    } else {
      stop("unrecognized generator specified");
    }
  }
  
  // step = s * (L * z)
  NumericVector step(d);
  if (d > 1) {
    for (int r = 0; r < d; ++r) {
      double acc = 0.0;
      for (int c = 0; c <= r; ++c) acc += this->L(r, c) * z[c];
      step[r] = this->s * this->D[r] * acc;
    }
  } else {
    step[0] = this->s * z[0];
  }
  
  for (int u = 0; u < d; ++u) w_new[u] = w_old_block[u] + step[u];
  
  Verbosity::instance().log(TRACE, 
                "[corrBlock] proposal w=", w_new);
  
  return w_new;
}

// adapt s and covariance
void CorrBlock::adapt(const double a_hat,
                      const std::vector<std::shared_ptr<StateParameter>>& parameters,
                      const int m) {
  
  // do not adapt once frozen
  if (this->frozen) return;
  
  // if iteration is beyond requested warmup period, freeze
  if (m >= this->warmup) {
    this->s = eps_adapter.freeze();
    this->cov_adapter.freeze();
    this->frozen = true;
    return;
  }
  
  // non-finite probabilities cannot be used in adaptation
  if (!std::isfinite(a_hat)) return;
  
  bool run_s;
  bool run_covD;
  bool run_covL;
  if (!this->learn_L && !this->learn_D) {
    run_s = true;
    run_covD = false;
    run_covL = false;
  } else if (this->learn_L) {
    if (m <= warmup / 5) {
      run_s = true;
      run_covD = false;
      run_covL = false;
    } else if (m <= warmup / 5 * 2) {
      run_s = false;
      run_covD = true;
      run_covL = false;
    } else if (m <= warmup / 5 * 4) {
      run_s = false;
      run_covD = false;
      run_covL = true;
    } else {
      run_s = true;
      run_covD = false;
      run_covL = false;
    }
  } else {
    if (m <= warmup / 3) {
      run_s = true;
      run_covD = false;
      run_covL = false;
    } else if (m <= warmup / 3 * 2) {
      run_s = false;
      run_covD = true;
      run_covL = false;
    } else {
      run_s = true;
      run_covD = false;
      run_covL = false;
    }
  }
  
  // adapt s
  if (this->learn_s && run_s) this->s = eps_adapter.update(a_hat);
  
  // if learning L or D, convert parameters to working/update scale
  // and pass to covariance adapter
  // warmup / 5 cut-off is to allow time for a "reasonable" value of s before updating
  // covariance
  if ((this->learn_D || this->learn_L) && (run_covD || run_covL) && (!cov_adapter.isFrozen())) {
    
    NumericVector w(this->d);
    for (int u = 0; u < this->d; ++u) {
      const double theta = parameters[u]->value;
      w[u] = parameters[u]->to_working(theta);
    }
    cov_adapter.update(w, this->L, this->D, m, run_covD, run_covL);
  }
  
  
  // keep progression information
  this->adapter_track(a_hat);
  
}
