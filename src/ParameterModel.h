#ifndef ParameterModel_H
#define ParameterModel_H

#include <Rcpp.h>
#include <memory>
#include "Pkg_types.h"

class StateParameter;

class ParameterModel {
public:
  
  // string description of model type
  Pkg::Model typeid_ = Pkg::Model::UNK;
  
  // string description of the distribution to which the model pertains
  std::string dist_param_name = "unk";
  
  // all parameters required of the model
  std::vector<std::shared_ptr<StateParameter>> parameters;
  
  // compute the distribution parameter at the provided t
  double (*fast_scalar)(const ParameterModel*, const double);
  Rcpp::NumericVector (*fast_vector)(const ParameterModel*, const Rcpp::NumericVector&);
  
  inline double par_from_tau(double t) const {
    return fast_scalar(this, t);
  }
  
  inline Rcpp::NumericVector par_from_tau_vec(const Rcpp::NumericVector& t) const {
    return fast_vector(this, t);
  }
  
  inline void par_from_tau(const Rcpp::NumericVector& tau, double* out) const {
    for (int i = 0; i < tau.size(); ++i)
      out[i] = fast_scalar(this, tau[i]);
  }
  
  // provided a vector of parameter values and t, evaluate the model
  // this is required for the Normalization calculated
  static double computeParameter(const Rcpp::NumericVector& params, double t) = delete;
  
  static Rcpp::NumericVector computeParameter(const Rcpp::NumericVector& params, 
                                              const Rcpp::NumericVector& t) = delete;
  
  // retrieve all current parameter values
  virtual Rcpp::NumericVector getParameters() const;
  
  // set all parameter values
  virtual void setParameters(const Rcpp::NumericVector& newParameterValues);
  
  ParameterModel(Rcpp::List model, std::string dist_param_name_, 
                 Pkg::Components component_,
                 Pkg::Model typeid_,
                 std::vector<std::shared_ptr<StateParameter>>& parameters);

  
  virtual std::unique_ptr<ParameterModel> clone() const = 0;
  virtual ~ParameterModel() = default;
  
  ParameterModel(ParameterModel&&) noexcept = default;
  ParameterModel& operator=(ParameterModel&&) noexcept = default;
  
protected:
  ParameterModel(const ParameterModel& other);
  ParameterModel() = default;
};


class ConstantParameterModel : public ParameterModel {
public:
  
  static double fast_fn_scalar(const ParameterModel* m, const double t);
  static Rcpp::NumericVector fast_fn_vector(const ParameterModel* m, const Rcpp::NumericVector& t);
  
  static double computeParameter(const Rcpp::NumericVector& params, double t);
  static Rcpp::NumericVector computeParameter(const Rcpp::NumericVector& params, 
                                              const Rcpp::NumericVector& t);
  
  ConstantParameterModel() : ParameterModel() {}
  
  ConstantParameterModel(Rcpp::List model, std::string dist_param_name_,
                         Pkg::Components component_,
                         std::vector<std::shared_ptr<StateParameter>>& parameters) :
    ParameterModel(model, dist_param_name_, component_, Pkg::Model::CONSTANT, parameters) {
    if (this->parameters.size() != 1) Rcpp::stop("wrong parameter count provided for constant model");
    fast_scalar = &ConstantParameterModel::fast_fn_scalar;
    fast_vector = &ConstantParameterModel::fast_fn_vector;
  };
  
  std::unique_ptr<ParameterModel> clone() const override {
    return std::make_unique<ConstantParameterModel>(*this);
  }
  
};

class SigmoidParameterModel : public ParameterModel {
public:

  static double fast_fn_scalar(const ParameterModel* m, const double t);
  static Rcpp::NumericVector fast_fn_vector(const ParameterModel* m, const Rcpp::NumericVector& t);
  
  static double computeParameter(const Rcpp::NumericVector& params, double t);
  static Rcpp::NumericVector computeParameter(const Rcpp::NumericVector& params, 
                                              const Rcpp::NumericVector& t);
  
  SigmoidParameterModel() : ParameterModel() {}
  
  SigmoidParameterModel(Rcpp::List model, std::string dist_param_name_,
                        Pkg::Components component_,
                        std::vector<std::shared_ptr<StateParameter>>& parameters) :
    ParameterModel(model, dist_param_name_, component_, Pkg::Model::SIGMOID, parameters) {
    if (this->parameters.size() != 4) Rcpp::stop("wrong parameter count provided for sigmoid model");
    fast_scalar = &SigmoidParameterModel::fast_fn_scalar;
    fast_vector = &SigmoidParameterModel::fast_fn_vector;
  };
  
  std::unique_ptr<ParameterModel> clone() const override {
    return std::make_unique<SigmoidParameterModel>(*this);
  }
};

class LinearParameterModel : public ParameterModel {
public:
  
  static double fast_fn_scalar(const ParameterModel* m, const double t);
  static Rcpp::NumericVector fast_fn_vector(const ParameterModel* m, const Rcpp::NumericVector& t);
  
  static double computeParameter(const Rcpp::NumericVector& params, double t);
  static Rcpp::NumericVector computeParameter(const Rcpp::NumericVector& params, 
                                              const Rcpp::NumericVector& t);
  
  LinearParameterModel() : ParameterModel() {}
  
  LinearParameterModel(Rcpp::List model, std::string dist_param_name_, 
                       Pkg::Components component_,
                       std::vector<std::shared_ptr<StateParameter>>& parameters) :
    ParameterModel(model, dist_param_name_, component_, Pkg::Model::LINEAR, parameters) {
    if (this->parameters.size() != 3) Rcpp::stop("wrong parameter count provided for linear model");
    fast_scalar = &LinearParameterModel::fast_fn_scalar;
    fast_vector = &LinearParameterModel::fast_fn_vector;
  };
  
  std::unique_ptr<ParameterModel> clone() const override {
    return std::make_unique<LinearParameterModel>(*this);
  }
};

class ExpParameterModel : public ParameterModel {
public:
  
  static double fast_fn_scalar(const ParameterModel* m, const double t);
  static Rcpp::NumericVector fast_fn_vector(const ParameterModel* m, const Rcpp::NumericVector& t);
  
  static double computeParameter(const Rcpp::NumericVector& params, double t);
  static Rcpp::NumericVector computeParameter(const Rcpp::NumericVector& params, 
                                              const Rcpp::NumericVector& t);
  
  ExpParameterModel() : ParameterModel() {}
  
  ExpParameterModel(Rcpp::List model, std::string dist_param_name_,
                    Pkg::Components component_,
                    std::vector<std::shared_ptr<StateParameter>>& parameters) :
    ParameterModel(model, dist_param_name_, component_, Pkg::Model::EXP, parameters) {
    if (this->parameters.size() != 3) Rcpp::stop("wrong parameter count provided for exp model");
    fast_scalar = &ExpParameterModel::fast_fn_scalar;
    fast_vector = &ExpParameterModel::fast_fn_vector;
  };
  
  std::unique_ptr<ParameterModel> clone() const override {
    return std::make_unique<ExpParameterModel>(*this);
  }
};

class StepParameterModel : public ParameterModel {
public:
  
  static double fast_fn_scalar(const ParameterModel* m, const double t);
  static Rcpp::NumericVector fast_fn_vector(const ParameterModel* m, const Rcpp::NumericVector& t);
  
  static double computeParameter(const Rcpp::NumericVector& params, double t);
  static Rcpp::NumericVector computeParameter(const Rcpp::NumericVector& params, 
                                              const Rcpp::NumericVector& t);
  
  StepParameterModel() : ParameterModel() {}
  
  StepParameterModel(Rcpp::List model,  std::string dist_param_name_,
                     Pkg::Components component_,
                     std::vector<std::shared_ptr<StateParameter>>& parameters) :
    ParameterModel(model, dist_param_name_, component_, Pkg::Model::STEP, parameters) {
    if (this->parameters.size() != 3) Rcpp::stop("wrong parameter count provided for step model");
    fast_scalar = &StepParameterModel::fast_fn_scalar;
    fast_vector = &StepParameterModel::fast_fn_vector;
  };
  
  std::unique_ptr<ParameterModel> clone() const override {
    return std::make_unique<StepParameterModel>(*this);
  }
};

#endif
