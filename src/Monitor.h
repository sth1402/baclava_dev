#ifndef MONITOR_H
#define MONITOR_H

#include <Rcpp.h>
#include "Compartment.h"
#include "Step_MH.h"
#include "Step_Gibbs.h"

class Monitor {
public:
  Rcpp::CharacterVector param_names;
  
  Rcpp::NumericMatrix posterior;
  Rcpp::LogicalVector kACCEPT;
  double acceptance_prob;
  int acceptance_sum;
  int n_iterations;
  
  
  Monitor() {};
  Monitor(const Step_MH& group, int M_thin);
  Monitor(const Step_Gibbs& group, int M_thin);
  
  void update(const Step_MH& group, int ikept, bool keep);
  void update(const Step_Gibbs& group, int ikept, bool keep);
  
  void update_acceptance_prob();
  
};

#endif
