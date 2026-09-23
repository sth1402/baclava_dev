#include <Rcpp.h>
#include <cmath>
#include <memory>
#include "Distribution.h"
#include "ParameterModel.h"
#include "StateParameter.h"
#include "NegativeRateError.h"

using namespace Rcpp;

/*
 *  Distribution$toList = function() {
 *    params <- lapply(self$parameters, function(x) x$toList())
 *    names(params) <- names(self$parameters)
 *    list("meta_for" = "Distribution", 
 *         "type" = self$type,
 *         "params" = params)
 *  }
 */
Distribution::Distribution(List distribution, Pkg::Dist type_id,
                           Pkg::Components component_,
                           std::vector<std::shared_ptr<StateParameter>>& parameters) {
  
  if (!distribution.containsElementNamed("meta_for"))
    Rcpp::stop("Distribution constructor: provided distribution list invalid");
  
  if (Rcpp::as<std::string>(distribution["meta_for"]) != "Distribution")
    Rcpp::stop("Distribution constructor: provided list is not of type Distribution");
  
  this->typeid_ = type_id;
  
  try {
    if (!distribution.containsElementNamed("params")) 
      stop("distribution does not contain element params");
    List distribution_params = distribution["params"];
  
    if (!distribution_params.containsElementNamed("shape")) 
      stop("distribution params does not contain element shape");
    List shape_list = distribution_params["shape"];
  
    // shape will be a result from Model$toList() applied to a const_m
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
    
    if (!shape_list.containsElementNamed("meta_for"))
      Rcpp::stop("Distribution constructor: provided shape list invalid");
    
    if (Rcpp::as<std::string>(shape_list["meta_for"]) != "Model")
      Rcpp::stop("Distribution constructor: provided shape list is not of type Model");
    
    if (!shape_list.containsElementNamed("type") ||
        !shape_list.containsElementNamed("params"))
        Rcpp::stop("shape_list is not properly formatted");
    
    std::string shape_model_type = as<std::string>(shape_list["type"]);
    if (shape_model_type != "constant") stop("only constant shape models are currently supported");
  
    List shape_model_params = shape_list["params"];
    if (shape_model_params.size() != 1) stop("too many parameters provided in shape");
    auto shape_param_list = Rcpp::as<Rcpp::List>(shape_model_params[0]);
    
    StateParameter tmp_shape = StateParameter(shape_param_list, component_);
    
    bool already_there = false;
    for (auto p : parameters) {
      if (p->tag == tmp_shape.tag) {
        already_there = true;
        this->shape = p;
        break;
      }
    }

    if (!already_there) {
      parameters.push_back(std::make_shared<StateParameter>(tmp_shape));
      this->shape = parameters.back();
    }    

    if (!distribution_params.containsElementNamed("mu") && 
        !distribution_params.containsElementNamed("rate") && 
        !distribution_params.containsElementNamed("scale") && 
        !distribution_params.containsElementNamed("median")) 
        stop("distribution parameters does not contain element mu, rate, median, or scale");
    
    int too_many = static_cast<int>(distribution_params.containsElementNamed("mu")) +
      static_cast<int>(distribution_params.containsElementNamed("rate")) +
      static_cast<int>(distribution_params.containsElementNamed("scale")) +
      static_cast<int>(distribution_params.containsElementNamed("median"));
    if (too_many > 1) stop("only 1 of mu/rate/scale/median can be provided");
    
    if (distribution_params.containsElementNamed("mu")) this->second_param = Pkg::Param::MU;
    if (distribution_params.containsElementNamed("rate")) this->second_param = Pkg::Param::RATE;
    if (distribution_params.containsElementNamed("scale")) this->second_param = Pkg::Param::SCALE;
    if (distribution_params.containsElementNamed("median")) this->second_param = Pkg::Param::MEDIAN;
    
  } catch (const std::exception& e) {
    Rcpp::stop("Distribution constructor: missing or malformed parameter list: %s", e.what());
  }
  
  // We do not store the other parameter here -- it is housed alongside
  // a Distribution in Compartment
}


double Distribution::getParameters() const {
  if (!shape) 
    Rcpp::stop("Distribution::getParameters: shape is null");
  return shape->value;
};
