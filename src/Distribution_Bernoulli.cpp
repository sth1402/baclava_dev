#include <Rcpp.h>
#include <cmath>
#include <memory>
#include "Distribution_Bernoulli.h"
#include "ParameterModel.h"
#include "StateParameter.h"
#include "NegativeRateError.h"

using namespace Rcpp;

/* Provided the ParameterModel defining the model of p and a current value
   of tau, returns the evaluated probability */
double Distribution_Bernoulli::get_prob(const double tau, 
                                       const ParameterModel* parameterModel) const {
  
  // evaluate distributional parameter under current functional parameters
  double individual_prob = parameterModel->par_from_tau(tau);
  
  if (individual_prob < 0.0 || individual_prob > 1.0) 
    throw NegativeRateError("Parameters yield an out of range probability");
  
  // clamp
  individual_prob = std::min(1.0 - 1e-8, std::max(1e-8, individual_prob));
  
  return individual_prob;
  
}

/* calculate density at x given current tau */
double Distribution_Bernoulli::dfuncf(
    const double x, const double tau, 
    bool log, 
    const ParameterModel* parameterModel)  const {
  
  // evaluate distributional parameter under current functional parameters
  double individual_prob = this->get_prob(tau, parameterModel);
  
  // evaluate density at x given current distributional parameters
  return R::dbinom(x, 1, individual_prob, log);
}

double Distribution_Bernoulli::fast_dfuncf_scalar(
    const Distribution* dist,
    const double x, const double tau, 
    bool log, 
    const ParameterModel* parameterModel) {
  return dist->dfuncf(x, tau, log, parameterModel);
}

NumericVector Distribution_Bernoulli::fast_dfuncf_vector(
    const Distribution* dist,
    const NumericVector& x, const NumericVector& tau, 
    bool log, 
    const ParameterModel* parameterModel) {
  NumericVector result(x.size());
  for(int i = 0; i < x.size(); ++i) 
    result[i] = dist->dfuncf(x[i], tau[i], log, parameterModel);
  return result;
}

/* calculate distribution function at x given current tau */
double Distribution_Bernoulli::pfuncf(
    const double x, const double tau, 
    bool lower, bool log, 
    const ParameterModel* parameterModel)  const {
  
  // evaluate distributional parameter under current functional parameters
  double individual_prob = this->get_prob(tau, parameterModel);
  
  return R::pbinom(x, 1, individual_prob, lower, log);
}

double Distribution_Bernoulli::fast_pfuncf_scalar(
    const Distribution* dist,
    const double x, const double tau, 
    bool lower, bool log, 
    const ParameterModel* parameterModel) {
  return dist->pfuncf(x, tau, lower, log, parameterModel);
}

NumericVector Distribution_Bernoulli::fast_pfuncf_vector(
    const Distribution* dist,
    const NumericVector& x, const NumericVector& tau, 
    bool lower, bool log, 
    const ParameterModel* parameterModel) {
  NumericVector result(x.size());
  for(int i = 0; i < x.size(); ++i) 
    result[i] = dist->pfuncf(x[i], tau[i], lower, log, parameterModel);
  return result;
}

/* calculate quantile function at x given current tau */
double Distribution_Bernoulli::qfuncf(
    const double x, const double tau, 
    bool lower, bool log, 
    const ParameterModel* parameterModel)  const {
  
  // evaluate distributional parameter under current functional parameters
  double individual_prob = this->get_prob(tau, parameterModel);
  
  return R::qbinom(x, 1, individual_prob, lower, log);
}
