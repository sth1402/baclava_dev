#include <Rcpp.h>
#include <cmath>
#include <memory>
#include "Distribution_Weibull.h"
#include "ParameterModel.h"
#include "StateParameter.h"
#include "NegativeRateError.h"

using namespace Rcpp;

/* constructor */
Distribution_Weibull::Distribution_Weibull(const Distribution_Weibull& other) : Distribution(other) {
  if (other.shape != nullptr) this->shape = other.shape;
}

/* 
 * convert functional parameters to distribution's scale;
 * param: current values of the functional parameters
*/
double Distribution_Weibull::parameter2scale(const double param) const {
  
  const double shape_val = shape ? shape->value : NA_REAL;
  
  if (scale_fn == nullptr) stop("scale_fn is not defined");
  return scale_fn(param, shape_val);
}

/* calculate density at x given current tau */
double Distribution_Weibull::dfuncf(
    const double x, const double tau, 
    bool log, 
    const ParameterModel* parameterModel)  const {
  if (x < 0.0) throw NegativeRateError("negative value provided to dweibull");
  
  // evaluate structural parameter under current functional parameters
  double structural_param = parameterModel->par_from_tau(tau);
  
  // convert evaluated structural parameter to distribution scale parameter
  double dist_scales = parameter2scale(structural_param);
  
  // evaluate density at x given current distributional parameters
  return R::dweibull(x, shape->value, dist_scales, log);
}

double Distribution_Weibull::fast_dfuncf_scalar(
    const Distribution* dist,
    const double x, const double tau, 
    bool log, 
    const ParameterModel* parameterModel) {
  return dist->dfuncf(x, tau, log, parameterModel);
}

NumericVector Distribution_Weibull::fast_dfuncf_vector(
    const Distribution* dist,
    const NumericVector& x, const NumericVector& tau, 
    bool log, 
    const ParameterModel* parameterModel) {
  NumericVector result(x.size());
  for(int i = 0; i < x.size(); ++i) result[i] = 
    dist->dfuncf(x[i], tau[i], log, parameterModel);
  return result;
}

/* calculate distribution function at x given current tau */
double Distribution_Weibull::pfuncf(
    const double x, const double tau, 
    bool lower, bool log, 
    const ParameterModel* parameterModel)  const {
  
  if (x < 0.0) throw NegativeRateError("negative value provided to pweibull");
  
  // evaluate structural parameter under current functional parameters
  double structural_param = parameterModel->par_from_tau(tau);
  
  // convert evaluated structural parameter to distribution scale parameter
  double dist_scales = parameter2scale(structural_param);

  // evaluate distribution at x given current distributional parameters
  return R::pweibull(x, shape->value, dist_scales, lower, log);
}

double Distribution_Weibull::fast_pfuncf_scalar(
    const Distribution* dist,
    const double x, const double tau, 
    bool lower, bool log, 
    const ParameterModel* parameterModel) {
  return dist->pfuncf(x, tau, lower, log, parameterModel);
}

NumericVector Distribution_Weibull::fast_pfuncf_vector(
    const Distribution* dist,
    const NumericVector& x, const NumericVector& tau, 
    bool lower, bool log, 
    const ParameterModel* parameterModel) {
  NumericVector result(x.size());
  for(int i = 0; i < x.size(); ++i) result[i] = 
    dist->pfuncf(x[i], tau[i], lower, log, parameterModel);
  return result;
}

/* calculate quantile function at x given current tau */
double Distribution_Weibull::qfuncf(
    const double x, const double tau, 
    bool lower, bool log, 
    const ParameterModel* parameterModel)  const {
  
  if (x < 0.0) throw NegativeRateError("negative value provided to qweibull");
  
  // evaluate structural parameter under current functional parameters
  double structural_param = parameterModel->par_from_tau(tau);
  
  // convert evaluated structural parameter to distribution scale parameter
  double dist_scales = parameter2scale(structural_param);
  
  // evaluate quantile function at x given current distributional parameters
  double result = R::qweibull(x, shape->value, dist_scales, lower, log);
  return result;
}

/*
 * Takes the age-dependent structural parameter, converts it to a standard 
 *   Weibull scale, and, when the shape indicates non-exponential behavior, 
 *   applies the analytic correction that produces the curvature-stabilized 
 *   "optimal" scale.
 */
double Distribution_Weibull::compute_optimal_scale(
    const double tau, 
    const ParameterModel* parameterModel) const {
  
  // evaluate structural parameter under current functional parameters
  double structural_param = parameterModel->par_from_tau(tau);
  
  // convert evaluated structural parameter to distribution scale parameter
  double dist_scales = parameter2scale(structural_param);
  
  // k ~ 1 => exponential; no adjustment in scale-space.
  if (shape->value < 1.0001) return dist_scales;
  
  // Convert scale -> rate (= 1/s)
  double rate_vec = 1.0 / dist_scales;
  
  // correct scale
  const double d = shape->value / (shape->value - 1.0);
  const double tmp = 1.0 / (1.0 - shape->value);
  
  double b1 = - std::pow(shape->value, shape->value * tmp) * std::pow(rate_vec, tmp);
  double b2 =   std::pow(shape->value, tmp) * std::pow(rate_vec, tmp);
  
  double lambda_opt = std::pow((b1 + b2) * d, -1.0 / d);
  double scale_opt  = 1.0 / lambda_opt;
  
  return scale_opt;
}
