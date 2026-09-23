#include <Rcpp.h>
#include <cmath>
#include "ParameterModel.h"
#include "StateParameter.h"
#include "NegativeRateError.h"

using namespace Rcpp;

namespace {
void verify_positive(const NumericVector& result) {
  for (double v : result) {
    if (v < 0.0 || !std::isfinite(v)) throw NegativeRateError("current model yields negative distribution parameter");
  }
}

void verify_positive(double result) {
  if (result < 0.0 || !std::isfinite(result)) throw NegativeRateError("current model yields negative distribution parameter");
}

}

ParameterModel::ParameterModel(const ParameterModel& other) {
  typeid_ = other.typeid_;
  dist_param_name = other.dist_param_name;
  
  parameters.reserve(other.parameters.size());
  for (const auto& p : other.parameters) {
    if (p) {
      parameters.push_back(std::make_shared<StateParameter>(*p));
    } else {
      parameters.push_back(nullptr);
    }
  }
}

/*
 * Model$toList = function() {
 *   param_list <- lapply(self$parameters, function(x) x$toList())
 *   param_tags <- lapply(self$parameters, "[[", "tag")
 *   names(param_list) <- param_tags
 *
 *   list("meta_for" = "Model",
 *        "type" = self$type, 
 *        "params" = param_list)
 * }
 */
ParameterModel::ParameterModel(List model, std::string dist_param_name_, 
                               Pkg::Components component_,
                               Pkg::Model typeid_,
                               std::vector<std::shared_ptr<StateParameter>>& parameters) {
  this->typeid_ = typeid_;
  this->dist_param_name = dist_param_name_;
  
  if (!model.containsElementNamed("meta_for")) 
    Rcpp::stop("ParameterModel constructor: provided parameter list invalid");

  if (Rcpp::as<std::string>(model["meta_for"]) != "Model")
    Rcpp::stop("ParameterModel constructor: provided list is not of type Model");
  
  if (!model.containsElementNamed("type") ||
      !model.containsElementNamed("params"))
      Rcpp::stop("model is not properly formatted");
  
  try {
    
    Rcpp::List list_params = as<Rcpp::List>(model["params"]);
    this->parameters.reserve(list_params.size());
    for (int i = 0; i < list_params.size(); i++) {
      StateParameter tmp_param = StateParameter(as<Rcpp::List>(list_params[i]), component_);
      
      bool already_there = false;
      for (auto p : parameters) {
        if (p->tag == tmp_param.tag) {
          already_there = true;
          this->parameters.push_back(p);
          break;
        }
      }
      
      if (!already_there) {
        parameters.push_back(std::make_shared<StateParameter>(tmp_param));
        this->parameters.push_back(parameters.back());
      }    
    }
    
  } catch (const std::exception& e) {
    Rcpp::stop("ParameterModel constructor: missing or malformed parameter list: %s", e.what());
  }
  
}


NumericVector ParameterModel::getParameters() const {
  const int k = parameters.size();
  NumericVector x(k);
  std::vector<std::string> nms;
  nms.reserve(k);
  
  for (int i = 0; i < k; ++i) {
    List tmp = parameters[i]->get_tag_value();
    x[i] = as<double>(tmp["value"]);
    nms.emplace_back(as<std::string>(tmp["tag"]));
  }
  
  std::transform(nms.begin(), nms.end(), nms.begin(),
                 [&](const std::string& s) { return dist_param_name + "." + s; });
  
  x.attr("names") = Rcpp::wrap(nms);
  return x;
}


// set all parameter values
void ParameterModel::setParameters(const Rcpp::NumericVector& newParameterValues) {
  if (newParameterValues.size() != static_cast<int>(parameters.size())) {
    Rcpp::stop("ParameterModel::setParameters: length mismatch");
  }
  
  for (int i = 0; i < parameters.size(); ++i) {
    if (parameters[i]) {
      parameters[i]->update(newParameterValues[i]);
    }
  }
}

// Constant 

double ConstantParameterModel::fast_fn_scalar(const ParameterModel* m, const double t) {
  return m->parameters[0]->value;
}

Rcpp::NumericVector ConstantParameterModel::fast_fn_vector(
    const ParameterModel* m,
    const Rcpp::NumericVector& t) {
  
  return Rcpp::NumericVector(t.size(), m->parameters[0]->value);
}

NumericVector ConstantParameterModel::computeParameter(const NumericVector& params,
                                                       const NumericVector& t) {
  return NumericVector(t.size(), params[0]);
}

double ConstantParameterModel::computeParameter(const NumericVector& params, double t) {
  return params[0];
}

// Sigmoid

NumericVector SigmoidParameterModel::fast_fn_vector(const ParameterModel* m,
                                                    const NumericVector& t) {
  if (std::fabs(m->parameters[3]->value) < 1e-12)
    throw NegativeRateError("current model divides by zero");
  NumericVector z = (t - m->parameters[2]->value) / m->parameters[3]->value;
  for (int i = 0; i < z.size(); ++i) {
    if (z[i] > 35.0) z[i] = 35.0;
    else if (z[i] < -35.0) z[i] = -35.0;
  }
  NumericVector result = m->parameters[0]->value + m->parameters[1]->value /
    (1.0 + exp(-z));
  verify_positive(result);
  return result;
}

double SigmoidParameterModel::fast_fn_scalar(const ParameterModel* m, const double t) {
  if (std::fabs(m->parameters[3]->value) < 1e-12)
    throw NegativeRateError("current model divides by zero");
  double z = (t - m->parameters[2]->value) / m->parameters[3]->value;
  z = std::max(std::min(z, 35.0), -35.0);
  double result = m->parameters[0]->value + m->parameters[1]->value /
    (1.0 + std::exp(-z));
  verify_positive(result);
  return result;
}

NumericVector SigmoidParameterModel::computeParameter(const NumericVector& params, 
                                                      const NumericVector& t) {
  if (params[3] < 1e-12 && params[3] > -1e-12)
    throw NegativeRateError("current model divides by zero");
  NumericVector z = (t - params[2]) / params[3];
  for (int i = 0; i < z.size(); ++i) {
    if (z[i] > 35.0) z[i] = 35.0;
    else if (z[i] < -35.0) z[i] = -35.0;
  }
  NumericVector result = params[0] + params[1] / (1.0 + exp(-z));
  verify_positive(result);
  return result;
}

double SigmoidParameterModel::computeParameter(const NumericVector& params, double t) {
  if (params[3] < 1e-12 && params[3] > -1e-12)
    throw NegativeRateError("current model divides by zero");
  double z = (t - params[2]) / params[3];
  z = std::max(std::min(z, 35.0), -35.0);
  double result = params[0] + params[1] / (1.0 + std::exp(-z));
       
  verify_positive(result);
  return result;
}

// Linear 

NumericVector LinearParameterModel::fast_fn_vector(const ParameterModel* m,const NumericVector& t) {
  NumericVector result = m->parameters[0]->value + m->parameters[1]->value * (t - m->parameters[2]->value);
  verify_positive(result);
  return result;
}

double LinearParameterModel::fast_fn_scalar(const ParameterModel* m, const double t) {
  double result = m->parameters[0]->value + m->parameters[1]->value * (t - m->parameters[2]->value);
  verify_positive(result);
  return result;
}

NumericVector LinearParameterModel::computeParameter(const NumericVector& params, 
                                                     const NumericVector& t) {
  NumericVector result = params[0] + params[1] * (t - params[2]);
  verify_positive(result);
  return result;
}

double LinearParameterModel::computeParameter(const NumericVector& params, double t) {
  double result = params[0] + params[1] * (t - params[2]);
  verify_positive(result);
  return result;
}

// Exponential Rate Model

NumericVector ExpParameterModel::fast_fn_vector(const ParameterModel* m, const NumericVector& t) {
  NumericVector z = m->parameters[1]->value * (t - m->parameters[2]->value);
  NumericVector result(t.size());
  for (int i = 0; i < z.size(); ++i) {
    if (z[i] > 700) throw NegativeRateError("exponential model overflow");  
    result[i] = m->parameters[0]->value * std::exp(z[i]);
  }
  verify_positive(result);
  return result;
}

double ExpParameterModel::fast_fn_scalar(const ParameterModel* m, const double t) {
  double z = m->parameters[1]->value * (t - m->parameters[2]->value);
  if (z > 700) throw NegativeRateError("exponential model overflow");
  double result = m->parameters[0]->value * std::exp(z);
  verify_positive(result);
  return result;
}

NumericVector ExpParameterModel::computeParameter(const NumericVector& params, 
                                                  const NumericVector& t) {
  NumericVector z = params[1] * (t - params[2]);
  NumericVector result(t.size());
  for (int i = 0; i < z.size(); ++i) {
    if (z[i] > 700) throw NegativeRateError("exponential model overflow");  
    result[i] = params[0] * std::exp(z[i]);
  }
  verify_positive(result);
  return result;
}

double ExpParameterModel::computeParameter(const NumericVector& params, double t) {
  double z = params[1] * (t - params[2]);
  if (z > 700) throw NegativeRateError("exponential model overflow");
  double result = params[0] * std::exp(z);
  verify_positive(result);
  return result;
}

// Step

NumericVector StepParameterModel::fast_fn_vector(const ParameterModel* m,const NumericVector& t) {
  NumericVector result(t.size());
  for (int i = 0; i < t.size(); ++i)
    result[i] = (t[i] < m->parameters[2]->value) ? m->parameters[0]->value : m->parameters[1]->value;
  verify_positive(result);
  return result;
}

double StepParameterModel::fast_fn_scalar(const ParameterModel* m, const double t) {
  double result = (t < m->parameters[2]->value) ?  m->parameters[0]->value : m->parameters[1]->value;
  verify_positive(result);
  return result;
}

NumericVector StepParameterModel::computeParameter(const NumericVector& params, 
                                                   const NumericVector& t) {
  NumericVector result(t.size());
  for (int i = 0; i < t.size(); ++i)
    result[i] = (t[i] < params[2]) ? params[0] : params[1];
  verify_positive(result);
  return result;
}

double StepParameterModel::computeParameter(const NumericVector& params, double t) {
  double result = (t < params[2]) ?  params[0] : params[1];
  verify_positive(result);
  return result;
}
