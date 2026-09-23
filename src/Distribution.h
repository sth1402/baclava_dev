#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <Rcpp.h>
#include <cmath>
#include <memory>
#include "Pkg_types.h"

class StateParameter;
class ParameterModel;

// Distribution functions including shape and density, distribution, and
// quantile functions; Also includes how to translate second parameter to
// scale as required by the specific distribution function type


class Distribution {
public:
  Pkg::Dist typeid_ = Pkg::Dist::UNK;
  
  // it is assumed that distribution will require a shape. 
  // allowed to be a NULL pointer if not required
  std::shared_ptr<StateParameter> shape;
  
  // specified the structural parameter; currently mu, rate, scale, or median
  Pkg::Param second_param = Pkg::Param::UNK;
  
  // constructors
  Distribution(Rcpp::List distribution, Pkg::Dist type_id, Pkg::Components component_,
               std::vector<std::shared_ptr<StateParameter>>& parameters);
  virtual ~Distribution() = default; 
  
  // density function
  virtual double dfuncf(const double x, const double tau, bool log, 
                        const ParameterModel* parameterModel) const = 0;
  
  virtual Rcpp::NumericVector dfuncf(const Rcpp::NumericVector& x, 
                                     const Rcpp::NumericVector& tau, 
                                     bool log, 
                                     const ParameterModel* parameterModel) const  {
    Rcpp::NumericVector result(x.size());
    for (int i = 0; i < x.size(); ++i) result[i] = this->dfuncf(x[i], tau[i], log, parameterModel);
    return result;
  };
  
  // distribution function
  virtual double pfuncf(const double x, const double tau, bool lower, bool log, 
                        const ParameterModel* parameterModel) const = 0;
  
  virtual Rcpp::NumericVector pfuncf(const Rcpp::NumericVector& x, 
                                     const Rcpp::NumericVector& tau, 
                                     bool lower, bool log, 
                                     const ParameterModel* parameterModel) const  {
    Rcpp::NumericVector result(x.size());
    for (int i = 0; i < x.size(); ++i) result[i] = this->pfuncf(x[i], tau[i], lower, log, parameterModel);
    return result;
  };
  
  // quantile function
  virtual double qfuncf(const double x, const double tau, bool lower, bool log, 
                        const ParameterModel* parameterModel) const = 0;
  
  virtual Rcpp::NumericVector qfuncf(const Rcpp::NumericVector& x, 
                                     const Rcpp::NumericVector& tau, 
                                     bool lower, bool log, 
                                     const ParameterModel* parameterModel) const  {
    Rcpp::NumericVector result(x.size());
    for (int i = 0; i < x.size(); ++i) result[i] = this->qfuncf(x[i], tau[i], lower, log, parameterModel);
    return result;
  };
  
  // retrieve distributional shape parameter
  virtual double getParameters() const;
  
  // convert structural parameter to distributional scale parameter using given
  // distributional shape and indicator of structural parameter type
  double (*scale_fn)(const double, const double);
  
  // convert structural parameter to distributional scale parameter
  virtual double parameter2scale(const double param) const = 0;
  virtual Rcpp::NumericVector parameter2scale(const Rcpp::NumericVector& param) const {
    Rcpp::NumericVector result(param.size());
    for (int i = 0; i < param.size(); ++i) result[i] = this->parameter2scale(param[i]);
    return result;
  };
  virtual void parameter2scale(double* param, double* out, int n) const {
    for (int i = 0; i < n; ++i) out[i] = this->parameter2scale(param[i]);
  };
  
  
  virtual std::unique_ptr<Distribution> clone() const = 0;

  // stabilize scale
  virtual double compute_optimal_scale(const double tau, 
                                       const ParameterModel* parameterModel) const = 0;
  
  virtual Rcpp::NumericVector compute_optimal_scale(const Rcpp::NumericVector& tau, 
                                                    const ParameterModel* parameterModel) const {
    Rcpp::NumericVector result(tau.size());
    for (int i = 0; i < tau.size(); ++i) result[i] = this->compute_optimal_scale(tau[i], parameterModel);
      return result;
  };
  
  // function to short-circuit to appropriate (Weibull/Gamma/Bernoulli) distributions
  double (*fast_dfunc_scalar)(const Distribution*, 
          const double x, const double tau, bool log, const ParameterModel*);
  Rcpp::NumericVector (*fast_dfunc_vector)(const Distribution*, 
                       const Rcpp::NumericVector& x, const Rcpp::NumericVector& tau, 
                       bool log, const ParameterModel*);
  
  double (*fast_pfunc_scalar)(const Distribution*, 
          const double x, const double tau, bool lower, bool log, const ParameterModel*);
  Rcpp::NumericVector (*fast_pfunc_vector)(const Distribution*, 
                       const Rcpp::NumericVector& x, const Rcpp::NumericVector& tau, 
          bool lower, bool log, const ParameterModel*);
  
  Distribution& operator=(const Distribution&) = delete;
  Distribution(Distribution&&) noexcept = default;
  Distribution& operator=(Distribution&&) noexcept = default;
  
protected:
  
  Distribution() = default;
  
  Distribution(const Distribution& other, int /*tag*/)
    : typeid_(other.typeid_),
      second_param(other.second_param),
      scale_fn(other.scale_fn){ }
  Distribution(const Distribution& other)
    : typeid_(other.typeid_),
      second_param(other.second_param),
      scale_fn(other.scale_fn){ }
  
};


#endif
