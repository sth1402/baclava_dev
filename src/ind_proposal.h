#ifndef IND_PROPOSAL_H
#define IND_PROPOSAL_H

#include <Rcpp.h>
#include "Pkg_types.h"

class Model;
struct EPGroup_State;

Rcpp::NumericVector compute_prob_indolent(const EPGroup_State& gs, 
                                          const Rcpp::NumericVector& tau,
                                          const Model* model);

Rcpp::IntegerVector rprop_indolent(const Pkg::EP endpoint_type, 
                                   const Rcpp::NumericVector& prob_indolent);

// vector of dimension equal to the number of participants
Rcpp::NumericVector dlog_prop_indolent(const Pkg::EP endpoint_type, 
                                       const Rcpp::NumericVector& prob_indolent,
                                       const Rcpp::IntegerVector& indolent);

#endif
