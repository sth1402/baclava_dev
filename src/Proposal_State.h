#ifndef PROPOSAL_STATE_H
#define PROPOSAL_STATE_H

#include <Rcpp.h>

/*
 *  convenience struct to compartmentalize current and new parameters,
 *  current parameters on the update (working) scale, propose new values, and 
 *  track auto rejection during an update step.
 */
class Proposal_State {
public:
  std::vector<double> parameters_old; // current parameters
  std::vector<double> parameters_new; // proposed parameters
  std::vector<double> w_old; // current parameters on update (working) scale
  bool auto_reject_block; // true = reject block always
  std::vector<bool> use_rnorm; // sampling function
  std::vector<bool> use_runif; 

  void init_from_group(const std::vector<std::shared_ptr<StateParameter>> parameters);
};

#endif
