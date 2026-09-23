#include "Monitor.h"
#include "Verbosity.h"
using namespace Rcpp;

Monitor::Monitor(const Step_MH& group, int M_thin) {
  
  n_iterations = 0;
  
  param_names = group.param_names;
  
  // trace matrices
  posterior = NumericMatrix(group.np, M_thin);
  posterior.attr("dimnames") = List::create(group.param_names, R_NilValue);

  kACCEPT = LogicalVector(M_thin);
  acceptance_prob = 0.0;
  acceptance_sum = 0;
}

Monitor::Monitor(const Step_Gibbs& group, int M_thin) {
  
  n_iterations = 0;
  
  param_names = group.param_names;
  
  // trace matrices
  posterior = NumericMatrix(group.np, M_thin);
  posterior.attr("dimnames") = List::create(group.param_names, R_NilValue);
  
  kACCEPT = LogicalVector(M_thin);
  acceptance_prob = 0.0;
  acceptance_sum = 0;
}

void Monitor::update(const Step_MH& group, int ikept, bool keep) {
  
  n_iterations++;
  
  Verbosity::instance().log(TRACE, "incrementing acceptance sum");

  // accumulate
  bool tmp_accepted = group.getAcceptStatus();
  if (tmp_accepted) acceptance_sum += 1;

  // store posterior and accept/reject only at the kept iterations
  // (past burnin and at thin)
  if (keep) {
    Verbosity::instance().log(TRACE, "saving posterior and accept/reject");
    if (ikept < posterior.ncol()) posterior(_, ikept) = group.getCurrentValues();
    if (ikept < kACCEPT.size()) kACCEPT(ikept) = group.getAcceptStatus();
  }
}

void Monitor::update(const Step_Gibbs& group, int ikept, bool keep) {
  
  n_iterations++;
  
  Verbosity::instance().log(TRACE, "incrementing acceptance sum");
  
  // accumulate
  bool tmp_accepted = group.getAcceptStatus();
  if (tmp_accepted) acceptance_sum += 1;
  
  // store posterior and accept/reject only at the kept iterations
  // (past burnin and at thin)
  if (keep) {
    Verbosity::instance().log(TRACE, "saving posterior and accept/reject");
    if (ikept < posterior.ncol()) posterior(_, ikept) = group.getCurrentValues();
    if (ikept < kACCEPT.size()) kACCEPT(ikept) = group.getAcceptStatus();
  }
}

void Monitor::update_acceptance_prob() {
  if (n_iterations <= 0) return;
  
  acceptance_prob = static_cast<double>(acceptance_sum) / 
      static_cast<double>(n_iterations);

}
