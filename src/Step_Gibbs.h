#ifndef STEP_GIBBS_H
#define STEP_GIBBS_H

#include <Rcpp.h>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include "Proposal_State.h"

class StateParameter;

class Step_Gibbs {
public:
  
  int np = 0;
  // vector to StateParameters updated in this MH step
  std::vector<std::shared_ptr<StateParameter>> parameters;
  
  // parameter names of the updated StateParameters
  Rcpp::CharacterVector param_names;
  
  
  Step_Gibbs(std::vector<std::shared_ptr<StateParameter>> params,
             Rcpp::CharacterVector param_names_)
    : np(static_cast<int>(params.size())),
      parameters(std::move(params)),
      param_names(std::move(param_names_)) {
    if (param_names.size() != static_cast<R_xlen_t>(parameters.size()))
      Rcpp::stop("Step_Gibbs: param_names size does not match parameters size");
  }
  
  Rcpp::NumericVector getCurrentValues() const;
  bool getAcceptStatus() const { return true; };
  
  Step_Gibbs() = default;
  Step_Gibbs(const Step_Gibbs&) = delete;
  Step_Gibbs& operator=(const Step_Gibbs&) = delete;
  
  Step_Gibbs(Step_Gibbs&& other) noexcept
    : np(other.np),
      parameters(std::move(other.parameters)),
      param_names(std::move(other.param_names)),
      params(std::move(other.params))
    {}
  
  Step_Gibbs& operator=(Step_Gibbs&& other) noexcept {
    if (this != &other) {
      np = other.np;
      parameters = std::move(other.parameters);
      param_names = std::move(other.param_names);
      params = std::move(other.params);
    }
    return *this;
  }
  
  ~Step_Gibbs() = default;
  
  Proposal_State params;
};

#endif
