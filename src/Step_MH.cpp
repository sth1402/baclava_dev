#include <Rcpp.h>
#include "StateParameter.h"
#include "Step_MH.h"
#include "CorrBlock.h"
#include "Verbosity.h"

using namespace Rcpp;

Step_MH::Step_MH(std::vector<std::shared_ptr<StateParameter>> params,
        Rcpp::CharacterVector param_names_)
  : np(static_cast<int>(params.size())),
    parameters(std::move(params)),
    param_names(std::move(param_names_)),
    accepted(false) {
  if (param_names.size() != static_cast<R_xlen_t>(parameters.size()))
    Rcpp::stop("Step_MH: param_names size does not match parameters size");
  
  /* binary flag indicating which model components are affected by proposal step*/
  impact = 0;
  for (auto p : parameters) impact |= static_cast<uint8_t>(p->component);
}


// draw a new update-scale proposal and map back to natural scale.
void Step_MH::propose() {
  
  // initialization of local struct -- sets vector lengths and initialzes
  // values
  this->proposal_state.init_from_group(this->parameters);
  
  // propose on working scale
  NumericVector w_old_block(this->np);
  for (int u = 0; u < np; ++u) w_old_block[u] = this->proposal_state.w_old[u];

  NumericVector w_new = this->corr_block.propose(w_old_block, this->proposal_state);

  // reflect on boundaries in working scale; map back to natural scale
  for (int u = 0; u < this->np; ++u) {
    
    double w = w_new[u];

    // map back to natural scale
    try {
      double theta_nat = this->parameters[u]->from_working(w);
      const double theta_new = this->parameters[u]->apply_boundaries(theta_nat);
      
      if (std::isnan(theta_new) || !R_finite(theta_new)) {
        throw std::runtime_error("NaN or non-finite theta_new");
      }
      this->proposal_state.parameters_new[u] = theta_new;
    } catch (...) {
      this->proposal_state.auto_reject_block = true;
    }

  }
  
  // if block auto flagged, revert all params in the block to old values
  if (this->proposal_state.auto_reject_block) {
    for (int u=0; u<np; ++u) 
      this->proposal_state.parameters_new[u] = this->proposal_state.parameters_old[u];
  }  
}  

/*
 * MH decision step
 * dlog_lik The difference in log-likelihood at the model level
 * m The current iteration
 * 
 * return a logical vector of accept/reject at the model level
 */

bool Step_MH::decide_commit_adapt(const double dlog_lik, int m)  {
  
  RNGScope scope;
  
  bool numerical_invalid = false;
  if (!std::isfinite(dlog_lik)) {
    Verbosity::instance().log(INFO, "NaN likelihood ");
    numerical_invalid = true;
  }

  // calculate priors at current and proposed parameter values
  NumericVector dlog_prior(this->np);
  for (int i = 0; i < this->np; ++i) {
    const double th_new = this->proposal_state.parameters_new[i];
    const double th_old = this->proposal_state.parameters_old[i];
    
    if (!this->parameters[i]) stop("Null shared_ptr at group->parameters[%d].", i);
    
    const double lp_new = this->parameters[i]->log_prior(th_new);
    const double lp_old = this->parameters[i]->log_prior(th_old);
    dlog_prior[i] = lp_new - lp_old;
    Verbosity::instance().log(TRACE, " lp_new=" + std::to_string(lp_new));
    Verbosity::instance().log(TRACE, " lp_old=" + std::to_string(lp_old));
  }

  // if parameter update scale does not align with scale of prior
  // add log |d theta / dw| term for transformed parameters (log/logit/logcv updates)
  NumericVector dlog_prop_jac(this->np);
  for (int i = 0; i < this->np; ++i) {
    const double th_new = this->proposal_state.parameters_new[i];
    const double th_old = this->proposal_state.parameters_old[i];
    const double logJ_new = this->parameters[i]->log_abs_dthetadw(th_new);
    const double logJ_old = this->parameters[i]->log_abs_dthetadw(th_old);
    Verbosity::instance().log(TRACE, " logJ_new=" + std::to_string(logJ_new));
    Verbosity::instance().log(TRACE, " logJ_old=" + std::to_string(logJ_old));
    dlog_prop_jac[i] = (logJ_new - logJ_old);
  }

  // priors are on the parameter level -- need them on the block level
  double prior_block = sum(dlog_prior);
  double jac_block = sum(dlog_prop_jac);

  for (int i = 0; i < this->np; ++i) {
    numerical_invalid = numerical_invalid || !std::isfinite(jac_block) ||
      !std::isfinite(prior_block) || !std::isfinite(dlog_lik);
  }
  
  double alpha_hat_block = 0.0;
  
  double MH_lr = dlog_lik + prior_block + jac_block;
  if (std::isnan(MH_lr)) MH_lr = R_NegInf;
  if (!std::isfinite(MH_lr)) MH_lr = R_NegInf;

  bool any_auto = this->proposal_state.auto_reject_block || numerical_invalid;

  const double logu = std::log(R::runif(0.0, 1.0));
  const bool accept_block = (!any_auto) && (logu < std::min(0.0, MH_lr));
    
  // parameter tracked by MC is the acceptance at the parameter level
  this->accepted = accept_block;
    
  // if the block is accepted, set all parameters to the proposed values
  if (accept_block) {
    for (int i = 0; i < this->np; ++i) {
      this->parameters[i]->update(this->proposal_state.parameters_new[i]);
    }
  } else {
    for (int i = 0; i < this->np; ++i) {
      this->parameters[i]->update(this->proposal_state.parameters_old[i]);
    }
  }

  // adaptation
  if (!this->corr_block.isFrozen() && !any_auto) {
    alpha_hat_block = std::exp(std::min(0.0, MH_lr));
    this->corr_block.adapt(alpha_hat_block, this->parameters, m);
  }
  
  return accept_block;
}

// retrieve current values of parameters being updated
NumericVector Step_MH::getCurrentValues() const {
  
  NumericVector value(this->np, 0.0);
  CharacterVector names(this->np);
  
  for (int i = 0; i < this->np; ++i) {
    value[i] = this->parameters[i]->value;
    names[i] = this->parameters[i]->tag;
  }
  value.names() = names;
  return value;
}

// set model parameters set to proposed values
void Step_MH::update_to_proposed() {
  for (int i = 0; i < np; ++i) this->parameters[i]->update(this->proposal_state.parameters_new[i]);
}

// set model parameters set to current values
void Step_MH::reset_to_current() {
  for (int i = 0; i < np; ++i) this->parameters[i]->update(this->proposal_state.parameters_old[i]);
}

// initialize proposal variables from list provided from R
void Step_MH::init_corr_block(const List& block) {
  this->corr_block = CorrBlock(block);
}

// extract proposal variables for return to R
List Step_MH::extract_epsilon() const {
  List out = this->corr_block.getFinalValues();
  out["param_names"] = this->param_names;
  return out;
}
