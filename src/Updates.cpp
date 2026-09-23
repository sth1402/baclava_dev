#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <algorithm>
#include "ind_proposal.h"
#include "tauHP_proposal.h"
#include "Step_MH.h"
#include "Model.h"
#include "NegativeRateError.h"
#include "Pkg_types.h"
#include "PriorFunctions.h"
#include "StateParameter.h"
#include "Updates.h"
#include "Verbosity.h"
#include "EPGroup_State.h"
#include "Step_Gibbs.h"

using namespace Rcpp;

/* General Metropolis Hastings Step */

void Updates::MH(std::vector<EPGroup_State>& groupstates,
                 Model* model, Step_MH& group, int m) {
  
  auto& vb = Verbosity::instance();
  if (vb.isEnabled(TRACE)) vb.log(TRACE, "cur ", group.getCurrentValues());
    
  // use cached values for current likelihood
  double loglik_cur = 0.0;
  for (auto& gs : groupstates) loglik_cur += sum(gs.log_likelihood());
    
  try {
    // propose new parameter values
    group.propose();
    
    // if auto-rejects (out of bounds, NA, etc, reject and keep current)
    if (group.proposal_state.auto_reject_block) {
      group.reset_to_current();
      group.setAcceptStatus(LogicalVector(group.np, false));
      return;
    }
    
    // update model parameters to the proposed values
    group.update_to_proposed();
      
    if (vb.isEnabled(TRACE)) vb.log(TRACE, "pro ", group.getCurrentValues());
      
    // likelihood under proposed model w.o. changing cache
    double loglik_new = 0.0;
    for (auto& gs : groupstates) loglik_new += sum(gs.log_likelihood(model, group.impact));
      
    if (vb.isEnabled(TRACE)) {
      vb.log(TRACE, "lik_cur ", loglik_cur);
      vb.log(TRACE, "lik_new ", loglik_new);
    }
      
    // accept/reject -- model is updated internally
    bool accept_block = group.decide_commit_adapt(loglik_new - loglik_cur, m);
      
    if (accept_block) {
      // if accepted, update impacted cached log-likelihoods
      for (auto& gs : groupstates) gs.update(model, group.impact);
    }

  } catch (const NegativeRateError& /*error*/) {
    group.reset_to_current();
    group.setAcceptStatus(LogicalVector(group.np, false));
  }
  
  if (vb.isEnabled(static_cast<VerboseLevel>(TRACE | INFO))) {
    vb.log(static_cast<VerboseLevel>(TRACE | INFO), "acc ", group.getCurrentValues());
  }
}

// joint update of psi and Z_indolent
void Updates::psiZ(std::vector<EPGroup_State>& groupstates, Model* model,
                   Step_MH& group, int m) {

  auto& vb = Verbosity::instance();
  if (vb.isEnabled(TRACE)) vb.log(TRACE, "cur ", group.getCurrentValues());

  // likelihood under current model and Z from cached values
  double loglik_cur = 0.0;
  for (auto& gs : groupstates) loglik_cur += sum(gs.log_likelihood());

  List prob_Z_cur(groupstates.size());
  // probability of indolence using current model
  for (int k = 0; k < groupstates.size(); ++k) {
    auto& gs = groupstates[k];
    if (gs.marginalize) continue;
    prob_Z_cur[k] = compute_prob_indolent(gs, gs.tau, model);
  }

  try {
    // propose new values
    group.propose();

    // if auto-rejects (out of bounds, NA, etc, reject and keep current)
    if (group.proposal_state.auto_reject_block) {
      group.reset_to_current();
      group.setAcceptStatus(LogicalVector(group.np, false));
      return;
    }

    // update model parameters to proposed values
    group.update_to_proposed();

    if (vb.isEnabled(TRACE)) vb.log(TRACE, "pro ", group.getCurrentValues());

    // likelihood under proposed model
    double loglik_new = 0.0;

    // proposal probabilities
    double prop_forward = 0.0;
    double prop_reverse = 0.0;

    // new latent indolence
    List Z_new(groupstates.size());
    
    for (int k = 0; k < groupstates.size(); ++k) {
      auto& gs = groupstates[k];
      
      // default impact does not include latent variable
      uint8_t aug_impact = group.impact;
      
      if (!gs.marginalize) {

        // latent probabilities and propose new values
        NumericVector prob_Z_new = compute_prob_indolent(gs, gs.tau, model);
        IntegerVector tmp_Z_new = rprop_indolent(gs.endpoint_type, prob_Z_new);

        // proposal probabilities
        // q(psi*, Z* | psi, Z) = q_psi(psi* | psi) * q_Z(Z* | psi*)
        // log q(Z_new | psi_new)
        prop_forward += sum(dlog_prop_indolent(gs.endpoint_type, prob_Z_new, tmp_Z_new));

        // q(psi, Z | psi*, Z*) = q_psi(psi | psi*) * q_Z(Z | psi)
        // log q(Z_cur | psi_cur)
        prop_reverse += sum(dlog_prop_indolent(gs.endpoint_type, prob_Z_cur[k], gs.Z));

        Z_new[k] = tmp_Z_new;

        // modify binary impact to include Z
        aug_impact |= static_cast<uint8_t>(Pkg::Components::Z);
        
        // likelihood calculated under proposed model with Z_new
        loglik_new += sum(gs.log_likelihood(model, tmp_Z_new, aug_impact));
      } else {
        Z_new[k] = gs.Z;
        // likelihood calculated under proposed model with Z_new
        loglik_new += sum(gs.log_likelihood(model, aug_impact));
      }

    }

    if (vb.isEnabled(TRACE)) {
      vb.log(TRACE, "psi_lik_cur ", loglik_cur);
      vb.log(TRACE, "psi_lik_new ", loglik_new);
      vb.log(TRACE, "psi_prop_for ", prop_forward);
      vb.log(TRACE, "psi_prop_rev ", prop_reverse);
    }

    // accept/reject -- model is updated internally
    bool accept_block =
      group.decide_commit_adapt(
          loglik_new - loglik_cur + (prop_reverse - prop_forward), m);

    if (accept_block) {
      for (int i = 0; i < groupstates.size(); ++i) {
        // if accepted, update cached log-likelihoods
        auto& gs = groupstates[i];
        
        uint8_t aug_impact = group.impact;
        if (gs.marginalize) {
          gs.update(model, aug_impact);
        } else {
          aug_impact |= static_cast<uint8_t>(Pkg::Components::Z);
          gs.update(model, as<IntegerVector>(Z_new[i]), aug_impact);
        }
      }
    }
  } catch (const NegativeRateError& /*error*/) {
    group.reset_to_current();
    group.setAcceptStatus(LogicalVector(group.np, false));
  }

  if (vb.isEnabled(static_cast<VerboseLevel>(TRACE | INFO))) {
    vb.log(static_cast<VerboseLevel>(TRACE | INFO), "acc ", group.getCurrentValues());
  }

}

// Gibbs update of beta
IntegerVector count_screens(const EPGroup_State& gs,
                            Model* model) {
  
  int n_screen_types = static_cast<int>(model->beta.size());
  IntegerVector missed(n_screen_types, 0);
  
  // count the total number of screens that occurred after the current
  // estimated age at time of healthy -> preclinical transition
  
  // no screens means that beta will not be updated
  if (static_cast<int>(gs.age_screen.size()) == 0) {
    return missed;
  }
  
  for (int i = 0; i < gs.age_screen.size(); ++i) {
    if (!std::isfinite(gs.tau[i])) continue;
    
    const auto& age_screen = gs.age_screen[i];
    int m = age_screen.size();
    if (gs.tau[i] > age_screen[m-1]) continue;
    
    const auto& screen_types = gs.screen_types[i];
    for (int j = 0; j < m; ++j) {
      if (age_screen[j] > gs.tau[i]) missed[screen_types[j]] += 1;
    }
  }
  return missed;
}

IntegerVector count_pos_screens(const EPGroup_State& gs,
                                Model* model) {
  
  int n_screen_types = static_cast<int>(model->beta.size());
  IntegerVector caught(n_screen_types, 0);
  
  if (gs.endpoint_type != Pkg::EP::PRECLINICAL) return caught;
  
  for (int i = 0; i < gs.screen_types.size(); ++i) {
    const auto& screen_types = gs.screen_types[i];
    int m = screen_types.size() - 1;
    int last_screen_type = screen_types[m];
    caught[last_screen_type] += 1;
  }
  return caught;
}

void Updates::beta(std::vector<EPGroup_State>& groupstates,
                   Model* model, Step_Gibbs& group) {
  
  auto& vb = Verbosity::instance();
  
  IntegerVector n_since_transition(static_cast<int>(model->beta.size()), 0);
  IntegerVector n_positive(static_cast<int>(model->beta.size()), 0);
  for (const auto& gs : groupstates) n_since_transition += count_screens(gs, model);
  for (const auto& gs : groupstates) n_positive += count_pos_screens(gs, model);
  
  if (vb.isEnabled(TRACE)) {
    vb.log(TRACE, "beta_cur ", group.getCurrentValues());
  }
  
  for (int i = 0; i < model->beta.size(); ++i) {
    if (model->beta[i].parameterModel->parameters[0]->fixed) continue;
    
    int success = n_positive[i];
    int total = n_since_transition[i];
    int failure = std::max(0, total - success);
    
    double new_beta = 0.0;
    const auto& prior = model->beta[i].parameterModel->parameters[0]->getPrior();
    
    if (prior.type() == "dbeta") {
      new_beta = R::rbeta(prior.params()[0] + success, prior.params()[1] + failure);
    } else if (prior.type() == "dunif") {
      new_beta = R::rbeta(1.0 + success, 1.0 + failure);
    } else {
      stop("unrecognized prior distribution");
    }
    
    model->beta[i].parameterModel->parameters[0]->value = new_beta;
  }
  
  // update cached log-likelihoods
  uint8_t aug_impact = static_cast<uint8_t>(Pkg::Components::BETA);
  for (auto& gs : groupstates) gs.update(model, aug_impact);

  if (vb.isEnabled(static_cast<VerboseLevel>(TRACE | INFO))) {
    vb.log(static_cast<VerboseLevel>(TRACE | INFO), "B_new ", group.getCurrentValues());
  }
}

void Updates::tau(std::vector<EPGroup_State>& groupstates, 
                  List& accepted, Model* model) {
  
  auto& vb = Verbosity::instance();
 
  uint8_t aug_impact = static_cast<uint8_t>(Pkg::Components::TAU);
 
  for (int i = 0; i < groupstates.size(); ++i) {
    auto& gs = groupstates[i];

    NumericVector endpoint_time = gs.endpoint_times;
    int n = gs.n;
    
    LogicalVector acc(n, false);
    accepted[i] = acc;
    
    // use cached likelihoods
    NumericVector loglik_cur = gs.log_likelihood();
    
    // store the current tau estimates
    NumericVector tau_cur = gs.tau;

    try {

      NumericVector wgt = compute_prob_tau_prep(gs, model);

      // compute interval-specific proposal probabilities under current model
      NumericVector prob_tau_cur = compute_prob_tau(gs, gs.tau, wgt, model);

      // propose new tau values using the interval-specific proposal probabilities
      NumericVector tau_new = rprop_tau(gs, gs.tau, prob_tau_cur, model);
      NumericVector prop_reverse = dlog_prop_tau(gs, tau_cur, prob_tau_cur, 
                                                 /*tau_ref=*/tau_new, model);
      
      // compute interval-specific proposal probabilities under proposed tau and model
      NumericVector prob_tau_new = compute_prob_tau(gs, tau_new, wgt, model);
      NumericVector prop_forward = dlog_prop_tau(gs, tau_new, prob_tau_new, 
                                                 /*tau_ref=*/tau_cur, model);
      
      // likelihood under proposed model returned per model w.o. touching cache
      NumericVector loglik_new = gs.log_likelihood(model, tau_new, aug_impact);
      
      if (vb.isEnabled(TRACE)) {
        vb.log(TRACE, "tau_lik_cur ", sum(loglik_cur));
        vb.log(TRACE, "tau_lik_new ", sum(loglik_new));
        vb.log(TRACE, "prop_revese ", sum(prop_reverse));
        vb.log(TRACE, "prop_forwrd ", sum(prop_forward));
      }
      
      NumericVector MH_logratio = loglik_new - loglik_cur + 
        (prop_reverse - prop_forward);
      MH_logratio[ is_na(MH_logratio) | is_nan(MH_logratio) ] = R_NegInf;

      for (int i = 0; i < n; ++i) {
        const double logu = std::log(R::runif(0.0, 1.0));
        acc[i] = (logu < std::min(0.0, MH_logratio[i]));
      }

      if (vb.isEnabled(static_cast<VerboseLevel>(TRACE | INFO))) {
        vb.log(static_cast<VerboseLevel>(TRACE | INFO), 
               "TMH " + std::to_string(Rcpp::min(MH_logratio)) + " "
                 + std::to_string(Rcpp::mean(MH_logratio)) + " "
                 + std::to_string(Rcpp::max(MH_logratio)) + " "
                 + std::to_string(sum(acc)));
      }
      
      if (is_true(any(acc))) {
        // Inf to Inf transitions are marked NA in acceptance vector
        LogicalVector tst = Rcpp::is_infinite(tau_new) & 
          Rcpp::is_infinite(tau_cur);
        
        // if any tau accepted, update cached log-likelihoods
        tau_cur[acc] = tau_new[acc];
        acc[tst] = NA_LOGICAL;
        
        gs.update(model, tau_cur, aug_impact);
      }

      if (vb.isEnabled(TRACE)) {
        vb.log(TRACE, "tau_lik_k ", sum(gs.log_likelihood()));
      }
      
    } catch (const NegativeRateError& error) {
      acc = LogicalVector(n, false);
    }

    accepted[i] = acc;

    if (vb.isEnabled(static_cast<VerboseLevel>(TRACE | INFO))) {
      vb.log(static_cast<VerboseLevel>(TRACE | INFO), 
             "T " + std::to_string(Rcpp::min(gs.tau)) + " "
               + std::to_string(Rcpp::mean(gs.tau)) + " "
               + std::to_string(Rcpp::max(gs.tau)));
    }
  }
  
}

