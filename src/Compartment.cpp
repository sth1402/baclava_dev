#include <Rcpp.h>
#include "Compartment.h"
#include "StateParameter.h"
#include "Distribution.h"
#include "EPGroup_State.h"
#include "Pkg_types.h"

using namespace Rcpp;


// retrieve all distribution parameters; assumes shape is last
NumericVector Compartment::getParameters() const {
  
  NumericVector dist_params = this->parameterModel->getParameters();
  
  NumericVector compartmentParameters = dist_params;
  if (distribution->shape) {
    CharacterVector names = compartmentParameters.names();
    compartmentParameters.push_back(this->distribution->shape->value);
    names.push_back(this->distribution->shape->tag);
    compartmentParameters.names() = names;
  }
  
  return compartmentParameters;
};

// reset all distribution parameters; assumes shape is last
void Compartment::setParameters(NumericVector newParameterValues) {
  
  if (this->distribution->shape) {
    this->distribution->shape->value = newParameterValues[newParameterValues.size() - 1];
  }
  
  this->parameterModel->setParameters(newParameterValues);
  
};

NumericVector Compartment::pfunc_ab(
    const NumericVector& a, const NumericVector& b, 
    const NumericVector& tau, bool return_log) const {
  
  if (this->distribution->typeid_ != Pkg::Dist::WEIBULL && 
      this->distribution->typeid_ != Pkg::Dist::GAMMA) {
    throw std::logic_error("pfunc_ab is only defined for continuous distributions (weibull/gamma).");
  }
  
  const ParameterModel* pm = parameterModel.get();
  const Distribution* dist = distribution.get();
  
  NumericVector logFA = dist->pfuncf(a, tau, 
                                     /*lower=*/true, /*log=*/true, 
                                     pm);
  NumericVector logFB = dist->pfuncf(b, tau, 
                                     /*lower=*/true, /*log=*/true, 
                                     pm);
      
  LogicalVector bad = (!Rcpp::is_finite(logFA) & !Rcpp::is_finite(logFB)) | (logFB <= logFA);
  NumericVector diff = logFA - logFB;
  diff[bad] = NA_REAL;
  NumericVector logInterval = logFB + Rcpp::log1p(-Rcpp::exp(diff));
  logInterval[bad | !Rcpp::is_finite(logInterval)] = R_NegInf;

  if (return_log) return logInterval;
  return exp(logInterval);
};

NumericVector Compartment::pfunc_ab_optimized(
    const EPGroup_State& gs, const NumericVector& tau, const double t0, 
    const Rcpp::NumericVector& wgt) const {
  
  int n = tau.size();
  std::vector<double> structural(n);
  std::vector<double> dist_scales(n);
  
  parameterModel->par_from_tau(tau, structural.data());
  distribution->parameter2scale(structural.data(), dist_scales.data(), n);
  
  const double shape_val = this->distribution->shape ? this->distribution->shape->value : NA_REAL;
  
  NumericVector result_probs = gs.pfunc_ab_optimized(dist_scales, 
                                                     this->distribution->typeid_, 
                                                     shape_val, t0,
                                                     wgt);

  return result_probs;
}

