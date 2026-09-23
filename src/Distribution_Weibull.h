#ifndef DISTRIBUTION_WEIBULL_H
#define DISTRIBUTION_WEIBULL_H

#include <Rcpp.h>
#include <cmath>
#include <memory>
#include <R_ext/Arith.h>
#include "Distribution.h"
#include "Pkg_types.h"

class StateParameter;
class ParameterModel;

class Distribution_Weibull : public Distribution {
public:
  
  // probability density function
  double dfuncf(const double x, const double tau, bool log, 
                const ParameterModel* parameterModel)  const override;
  
  // probability distribution function
  double pfuncf(const double x, const double tau, bool lower, bool log, 
                const ParameterModel* parameterModel)  const override;
  
  // probability quantile function
  double qfuncf(const double x, const double tau, bool lower, bool log, 
                const ParameterModel* parameterModel)  const override;
  
  // scale stabilization
  double compute_optimal_scale(const double tau, 
                               const ParameterModel* parameterModel) const override;
  
  // convert structural parameters mu/rate/scale/median to distribution's scale 
  // parameter given structural parameter (param), distribution's shape (shape_)
  
  // static inline double p2scale(const double param, const double shape_, const double t_ref, const double t0) {
  //   if (param <= 0.0 || param >= 1.0) return R_PosInf;
  //   if (t_ref <= t0) return R_PosInf;
  //   return (t_ref - t0) / std::pow(-std::log1p(-param), 1.0 / shape_);
  // }

  static inline double mu2scale(const double param, const double shape_) {
    double gamma_term = std::tgamma(1.0 + 1.0 / shape_);
    return param / gamma_term;
  }
  
  static inline double rate2scale(const double param, const double shape_) {
    return pow(param, -1.0 / shape_);
  }
  
  static inline double scale2scale(const double param, const double shape_) {
    return param;
  }
  
  static inline double median2scale(const double param, const double shape_) {
    return param / std::pow(std::log(2.0), 1.0 / shape_);
  }
  
  // combine provided structural parameter mu/rate/scale with distributional
  // parameter shape to obtain distributional scale parameter
  double parameter2scale(const double param) const override;
  
  Distribution_Weibull() : Distribution() {}
  
  Distribution_Weibull(Rcpp::List distribution, Pkg::Components component_,
                       std::vector<std::shared_ptr<StateParameter>>& parameters) : 
    Distribution(distribution, Pkg::Dist::WEIBULL, component_, parameters) {
    switch (second_param) {
    case Pkg::Param::MU:     scale_fn = &Distribution_Weibull::mu2scale; break;
    case Pkg::Param::RATE:   scale_fn = &Distribution_Weibull::rate2scale; break;
    case Pkg::Param::SCALE:  scale_fn = &Distribution_Weibull::scale2scale; break;
    case Pkg::Param::MEDIAN: scale_fn = &Distribution_Weibull::median2scale; break;
    default:     scale_fn = nullptr; break;
    }
    
    fast_dfunc_scalar = &Distribution_Weibull::fast_dfuncf_scalar;
    fast_dfunc_vector = &Distribution_Weibull::fast_dfuncf_vector;
    
    fast_pfunc_scalar = &Distribution_Weibull::fast_pfuncf_scalar;
    fast_pfunc_vector = &Distribution_Weibull::fast_pfuncf_vector;
  };
  
  Distribution_Weibull(const Distribution_Weibull& other);
  
  std::unique_ptr<Distribution> clone() const override {
    return std::make_unique<Distribution_Weibull>(*this);
  }
  
  static double fast_dfuncf_scalar(const Distribution* dist,
                                   const double x, const double tau, 
                                   bool log, 
                                   const ParameterModel* parameterModel);
  static double fast_pfuncf_scalar(const Distribution* dist,
                                   const double x, const double tau, 
                                   bool lower, bool log, 
                                   const ParameterModel* parameterModel);
  static Rcpp::NumericVector fast_dfuncf_vector(const Distribution* dist,
                                                const Rcpp::NumericVector& x, const Rcpp::NumericVector& tau, 
                                                bool log, 
                                                const ParameterModel* parameterModel);
  static Rcpp::NumericVector fast_pfuncf_vector(const Distribution* dist,
                                                const Rcpp::NumericVector& x, const Rcpp::NumericVector& tau, 
                                                bool lower, bool log, 
                                                const ParameterModel* parameterModel);
  
};

#endif
