#ifndef STEP_MH_H
#define STEP_MH_H

#include <Rcpp.h>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include "CorrBlock.h"
#include "Proposal_State.h"
#include "Pkg_types.h"

class StateParameter;

class Step_MH {
public:
  
  int np = 0;
  
  // vector of StateParameters updated in this MH step
  std::vector<std::shared_ptr<StateParameter>> parameters;
  
  // vector of model components this update impacts
  uint8_t impact = 0;
  
  // parameter names of the updated StateParameters
  Rcpp::CharacterVector param_names;

  Step_MH(std::vector<std::shared_ptr<StateParameter>> params,
          Rcpp::CharacterVector param_names_);
  
  void propose();
  
  bool decide_commit_adapt(const double dlog_lik, int m);
  
  void update_to_proposed();
  void reset_to_current();
  
  Rcpp::NumericVector getCurrentValues() const;
  bool getAcceptStatus() const { return accepted; };
  void setAcceptStatus(bool accepted_) { accepted = accepted_; }
  Rcpp::List getEvolution() const {return corr_block.getAdaptivePhase();}

  void init_corr_block(const Rcpp::List& block);
  
  Rcpp::List extract_epsilon() const;
  
  Step_MH() = default;
  Step_MH(const Step_MH&) = delete;
  Step_MH& operator=(const Step_MH&) = delete;
  
  Step_MH(Step_MH&& other) noexcept
    : np(other.np),
      parameters(std::move(other.parameters)),
      impact(std::move(other.impact)),
      param_names(std::move(other.param_names)),
      proposal_state(std::move(other.proposal_state)),
      accepted(std::move(other.accepted)),
      corr_block(std::move(other.corr_block))
  {}
  
  Step_MH& operator=(Step_MH&& other) noexcept {
    if (this != &other) {
      np = other.np;
      parameters = std::move(other.parameters);
      impact = std::move(other.impact);
      param_names = std::move(other.param_names);
      proposal_state = std::move(other.proposal_state);
      accepted = std::move(other.accepted);
      corr_block = std::move(other.corr_block);
    }
    return *this;
  }
  
  ~Step_MH() = default;
  
  Proposal_State proposal_state;
  
private:
  
  // indicating if each updated StateParameter is accepted/rejected
  bool accepted;

  // updating and adaptive procedures
  CorrBlock corr_block;
  
  static Rcpp::NumericMatrix diag(int d) {
    Rcpp::NumericMatrix I(d, d);
    for (int i = 0; i < d; ++i) I(i, i) = 1.0;
    return I;
  }
  
};

#endif
