#ifndef DISTRIBUTION_BERNOULLI_H
#define DISTRIBUTION_BERNOULLI_H

#include <Rcpp.h>
#include <cmath>
#include <memory>
#include "Distribution.h"
#include "Pkg_types.h"

class StateParameter;
class ParameterModel;

// Distribution functions including shape and density, distribution, and
// quantile functions; Also include how to translate second parameter to
// scale as required by the specific distribution function type

class Distribution_Bernoulli : public Distribution {
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
  
  // should never be used
  double compute_optimal_scale(
      const double tau, 
      const ParameterModel* parameterModel) const override {
        throw std::logic_error("Function not defined for this distribution type");
  };
  
  // should never be used
  static inline double mu2scale(const double param, const double shape_) {
    throw std::logic_error("Function not defined for this distribution type");;
  }
  static inline double rate2scale(const double param, const double shape_) {
    throw std::logic_error("Function not defined for this distribution type");;
  }
  static inline double scale2scale(const double param, const double shape_) {
    throw std::logic_error("Function not defined for this distribution type");;
  }
  static inline double median2scale(const double param, const double shape_) {
    throw std::logic_error("Function not defined for this distribution type");;
  }
  
  // should never be used
  double parameter2scale(const double param) const override {
    throw std::logic_error("Function not defined for this distribution type");;
  };
  
  // should never be used
  double getParameters() const override {
    throw std::logic_error("Function not defined for this distribution type");;
  };
  
  Distribution_Bernoulli() : Distribution() {}
  
  // has no shape parameter
  Distribution_Bernoulli(Rcpp::List distribution,
                         std::vector<std::shared_ptr<StateParameter>>& parameters) {
    this->typeid_ = Pkg::Dist::BERNOULLI;
    this->shape = nullptr;
    this->scale_fn = nullptr;
    this->second_param = Pkg::Param::UNK;
  };
  
  Distribution_Bernoulli(const Distribution_Bernoulli& other) : Distribution(other) {
    this->shape = nullptr;
    this->scale_fn = nullptr;
    this->second_param = Pkg::Param::UNK;
    this->typeid_ = other.typeid_;
    
    fast_dfunc_scalar = &Distribution_Bernoulli::fast_dfuncf_scalar;
    fast_dfunc_vector = &Distribution_Bernoulli::fast_dfuncf_vector;
    
    fast_pfunc_scalar = &Distribution_Bernoulli::fast_pfuncf_scalar;
    fast_pfunc_vector = &Distribution_Bernoulli::fast_pfuncf_vector;
  }
  
  std::unique_ptr<Distribution> clone() const override {
    return std::make_unique<Distribution_Bernoulli>(*this);
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
  
  double get_prob(const double tau, 
                  const ParameterModel* parameterModel) const;
  
};

#endif
