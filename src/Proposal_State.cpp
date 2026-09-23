#include <Rcpp.h>
#include "StateParameter.h"
#include "Proposal_State.h"

using namespace Rcpp;

void Proposal_State::init_from_group(const std::vector<std::shared_ptr<StateParameter>> parameters) {
  
  int p = static_cast<int>(parameters.size());
  
  // size and default
  parameters_old.resize(p);
  parameters_new.resize(p);
  w_old.resize(p);
  auto_reject_block = false;
  use_rnorm.resize(p);
  use_runif.resize(p);
  
  // set new and old as current value; convert old to update scale
  for (int i = 0; i < p; ++i) {
    const double theta = parameters[i]->get_value();
    parameters_old[i] = theta;
    parameters_new[i] = theta;
    w_old[i] = parameters[i]->to_working(theta);
    use_rnorm[i] = parameters[i]->prop_func == "rnorm";
    use_runif[i] = parameters[i]->prop_func == "runif";
  }
}
