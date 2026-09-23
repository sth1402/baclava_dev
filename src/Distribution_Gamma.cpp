#include <Rcpp.h>
#include <cmath>
#include <memory>
#include "Distribution_Gamma.h"
#include "ParameterModel.h"
#include "StateParameter.h"
#include "NegativeRateError.h"

using namespace Rcpp;

/* constructor */
Distribution_Gamma::Distribution_Gamma(const Distribution_Gamma& other) : Distribution(other) {
  if (other.shape != nullptr) this->shape = other.shape;//std::make_unique<StateParameter>(*other.shape);
};

/* 
 * given the structural parameter, combine with distributional parameter shape
 * to obtain the distributional scale parameter
 */
double Distribution_Gamma::parameter2scale(const double param) const {
  
  const double shape_val = shape ? shape->value : NA_REAL;
  
  if (!(shape_val > 0.0) || !R_finite(shape_val)) {
    return NA_REAL;
  }
  
  if (scale_fn == nullptr) stop("scale_fn is not defined");
  double result = scale_fn(param, shape_val);
  
  return result;
}

double Distribution_Gamma::get_scale(const double tau, 
                                     const ParameterModel* parameterModel) const {
  // evaluate structural parameter under current functional parameters
  double dist_parameter = parameterModel->par_from_tau(tau);
  
  // convert evaluated structural parameter to distribution scale parameter
  return this->parameter2scale(dist_parameter);
  
}

/* calculate density at x given current tau */
double Distribution_Gamma::dfuncf(
    const double x, const double tau, 
    bool log, 
    const ParameterModel* parameterModel)  const {
  
  if (x < 0.0) throw NegativeRateError("negative value provided to dweibull");

  // convert evaluated structural parameter to distribution scale parameter
  double dist_scales = get_scale(tau, parameterModel);
  
  // evaluate density at x given current distributional parameters
  return R::dgamma(x, shape->value, dist_scales, log);
}

double Distribution_Gamma::fast_dfuncf_scalar(
    const Distribution* dist,
    const double x, const double tau, 
    bool log, 
    const ParameterModel* parameterModel) {
  return dist->dfuncf(x, tau, log, parameterModel);
}

NumericVector Distribution_Gamma::fast_dfuncf_vector(
    const Distribution* dist,
    const NumericVector& x, const NumericVector& tau, 
    bool log, 
    const ParameterModel* parameterModel) {
  NumericVector result(x.size());
  for(int i = 0; i < x.size(); ++i) result[i] = 
    dist->dfuncf(x[i], tau[i], log, parameterModel);
  return result;
}

double Distribution_Gamma::pfuncf(
    const double x, const double tau, 
    bool lower, bool log, 
    const ParameterModel* parameterModel)  const {

  if (x < 0.0) throw NegativeRateError("negative value provided to dweibull");
  
  double dist_scales = get_scale(tau, parameterModel);
  
  return R::pgamma(x, shape->value, dist_scales, lower, log);
}
double Distribution_Gamma::fast_pfuncf_scalar(
    const Distribution* dist,
    const double x, const double tau, 
    bool lower, bool log, 
    const ParameterModel* parameterModel) {
  return dist->pfuncf(x, tau, lower, log, parameterModel);
}

NumericVector Distribution_Gamma::fast_pfuncf_vector(
    const Distribution* dist,
    const NumericVector& x, const NumericVector& tau, 
    bool lower, bool log, 
    const ParameterModel* parameterModel) {
  NumericVector result(x.size());
  for(int i = 0; i < x.size(); ++i) result[i] = 
    dist->pfuncf(x[i], tau[i], lower, log, parameterModel);
  return result;
}

double Distribution_Gamma::qfuncf(
    const double x, const double tau, 
    bool lower, bool log, 
    const ParameterModel* parameterModel)  const {
  
  if (x < 0.0) throw NegativeRateError("negative value provided to dweibull");

  double dist_scales = get_scale(tau, parameterModel);
  
  double result = R::qgamma(x, shape->value, dist_scales, lower, log);
  return result;
}

/*
 * Takes the age-dependent structural parameter, converts it to a standard 
 *   Gamma scale.
 */
double Distribution_Gamma::compute_optimal_scale(
    const double tau, 
    const ParameterModel* parameterModel) const {
  
  return get_scale(tau, parameterModel);
}
