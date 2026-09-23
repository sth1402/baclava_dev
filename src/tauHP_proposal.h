#ifndef TAUHP_PROPOSAL_H
#define TAUHP_PROPOSAL_H

#include <Rcpp.h>

class Model;
struct EPGroup_State;

Rcpp::NumericVector compute_prob_tau_prep(const EPGroup_State& gs, const Model* model);

Rcpp::NumericVector compute_prob_tau(const EPGroup_State& gs, 
                                     const Rcpp::NumericVector& tau, 
                                     const Rcpp::NumericVector& wgt, 
                                     const Model* model);

Rcpp::NumericVector rprop_tau(const EPGroup_State& gs, 
                              const Rcpp::NumericVector& tau, 
                              const Rcpp::NumericVector& prob_tau, 
                              const Model* model);

Rcpp::NumericVector dlog_prop_tau(const EPGroup_State& gs,
                                  const Rcpp::NumericVector& tau_eval,
                                  const Rcpp::NumericVector& prob_tau_eval,
                                  const Rcpp::NumericVector& tau_ref,
                                  const Model* model);


#endif
