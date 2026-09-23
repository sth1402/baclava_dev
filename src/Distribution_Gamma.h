#ifndef DISTRIBUTION_GAMMA_H
#define DISTRIBUTION_GAMMA_H

#include <Rcpp.h>
#include <cmath>
#include <memory>
#include "Distribution.h"
#include "Pkg_types.h"

class StateParameter;
class ParameterModel;

class Distribution_Gamma : public Distribution {
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
  
  // convert structural parameters mu/rate/scale to distribution scale 
  // parameter given structural parameter (param), distribution shape (shape_)
  // and a string indicating what the structural parameter is
  static inline double mu2scale(const double param, const double shape_) {
     return param / shape_;
  }
  
  static inline double rate2scale(const double param, const double shape_) {
    return 1.0 / param;
  }
  
  static inline double scale2scale(const double param, const double shape_) {
    return param;
  }
  
  static inline double median2scale(const double param, const double shape_) {
    return NA_REAL;
  }
  
  // combine provided structural parameter mu/rate/scale with distributional
  // parameter shape to obtain distributional scale parameter
  double parameter2scale(const double param) const override;
  
  Distribution_Gamma() : Distribution() {}
  
  Distribution_Gamma(Rcpp::List distribution, Pkg::Components component_,
                     std::vector<std::shared_ptr<StateParameter>>& parameters) : 
    Distribution(distribution, Pkg::Dist::GAMMA, component_, parameters) {
    switch (second_param) {
    case Pkg::Param::MU:     scale_fn = &Distribution_Gamma::mu2scale; break;
    case Pkg::Param::RATE:   scale_fn = &Distribution_Gamma::rate2scale; break;
    case Pkg::Param::SCALE:  scale_fn = &Distribution_Gamma::scale2scale; break;
    default:     scale_fn = nullptr; break;
    }
    
    fast_dfunc_scalar = &Distribution_Gamma::fast_dfuncf_scalar;
    fast_dfunc_vector = &Distribution_Gamma::fast_dfuncf_vector;

    fast_pfunc_scalar = &Distribution_Gamma::fast_pfuncf_scalar;
    fast_pfunc_vector = &Distribution_Gamma::fast_pfuncf_vector;
  };
  
  Distribution_Gamma(const Distribution_Gamma& other);
  
  std::unique_ptr<Distribution> clone() const override {
    return std::make_unique<Distribution_Gamma>(*this);
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
  
  
private:
  double get_scale(const double tau, 
                   const ParameterModel* parameterModel) const;
};

#endif
