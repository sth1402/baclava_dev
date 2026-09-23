#ifndef PARAMETERMODELFACTORY_H
#define PARAMETERMODELFACTORY_H
#include "ParameterModel.h"
#include "Pkg_types.h"

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
inline std::shared_ptr<ParameterModel>
  make_param_model_from_list(Pkg::Components component_,
                             const Rcpp::List& modeled_param,
                             const std::string& modeled_param_name,
                             std::vector<std::shared_ptr<StateParameter>>& parameters) {
    
    if (!modeled_param.containsElementNamed("meta_for"))
      Rcpp::stop("ParameterModel Factory: provided parameter list invalid");
    
    if (Rcpp::as<std::string>(modeled_param["meta_for"]) != "Model")
      Rcpp::stop("ParameterModel Factory: provided list is not of type Model");
    
    if (!modeled_param.containsElementNamed("type") ||
        !modeled_param.containsElementNamed("params"))
        Rcpp::stop("modeled_param is not properly formatted");
    
    const auto model_type = Rcpp::as<std::string>(modeled_param["type"]);

    if (model_type == "constant")  
      return std::make_shared<ConstantParameterModel>(modeled_param, 
                                                      modeled_param_name,
                                                      component_,
                                                      parameters);
    if (model_type == "sigmoid")   
      return std::make_shared<SigmoidParameterModel>(modeled_param, 
                                                     modeled_param_name,
                                                     component_,
                                                     parameters);
    if (model_type == "linear")    
      return std::make_shared<LinearParameterModel>(modeled_param, 
                                                    modeled_param_name,
                                                    component_,
                                                    parameters);
    if (model_type == "exp")       
      return std::make_shared<ExpParameterModel>(modeled_param, 
                                                 modeled_param_name,
                                                 component_,
                                                 parameters);
    if (model_type == "step")      
      return std::make_shared<StepParameterModel>(modeled_param, 
                                                  modeled_param_name,
                                                  component_,
                                                  parameters);
    
    throw std::invalid_argument("Parameter model type not recognized: " + model_type);
  }

/*
 *  Distribution$toList = function() {
 *    params <- lapply(self$parameters, function(x) x$toList())
 *    names(params) <- names(self$parameters)
 *    list("meta_for" = "Distribution", 
 *         "type" = self$type,
 *         "params" = params)
 *  }
 */
inline std::pair<Rcpp::List, std::string> extract_modeled_param(
    const Rcpp::List& d) {
  
  if (!d.containsElementNamed("meta_for"))
    Rcpp::stop("extract_modeled_param: provided distribution list invalid");
  
  if (Rcpp::as<std::string>(d["meta_for"]) != "Distribution")
    Rcpp::stop("extract_modeled_param: provided list is not of type Distribution");
  
  if (!d.containsElementNamed("params"))
      Rcpp::stop("Distribution missing 'params'");
    const Rcpp::List distribution_params = d["params"];
    
    // Only 1 parameter may be modeled
    if (distribution_params.containsElementNamed("mu"))
      return { distribution_params["mu"],   "mu"   };
    if (distribution_params.containsElementNamed("rate"))
      return { distribution_params["rate"], "rate" };
    if (distribution_params.containsElementNamed("scale"))
      return { distribution_params["scale"], "scale" };
    if (distribution_params.containsElementNamed("median"))
      return { distribution_params["median"], "median" };
    if (distribution_params.containsElementNamed("p"))
      return { distribution_params["p"],    "p"    };
    
    Rcpp::stop("Distribution parameters must include either 'mu', 'rate', 'median', or 'scale' (Weibull/Gamma) or 'p' (Bernoulli).");
  }
#endif
