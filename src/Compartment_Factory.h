#ifndef COMPARTMENT_FACTORY_H
#define COMPARTMENT_FACTORY_H

#include "Compartment.h"
#include "Distribution_Factory.h"
#include "ParameterModel_Factory.h"
#include "Verbosity.h"
#include "Pkg_types.h"

// Convenience functions to create Compartment objects

inline std::vector<Compartment> build_compartment(
    const Rcpp::List& dists, Pkg::Components component_,
    std::vector<std::shared_ptr<StateParameter>>& parameters) {
  
  std::vector<Compartment> cds;
  
  cds.reserve(dists.size());
  for (int i = 0; i < dists.size(); ++i) {
    
    const Rcpp::List dist_i = Rcpp::as<Rcpp::List>(dists[i]);
    /* dist_i
     *  Distribution$toList = function() {
     *    params <- lapply(self$parameters, function(x) x$toList())
     *    names(params) <- names(self$parameters)
     *    list("meta_for" = "Distribution", 
     *         "type" = self$type,
     *         "params" = params)
     *  }
     */
    
    if (!dist_i.containsElementNamed("meta_for"))
      Rcpp::stop("build_compartment: provided distribution list invalid");
    
    if (Rcpp::as<std::string>(dist_i["meta_for"]) != "Distribution")
      Rcpp::stop("build_compartment: provided list is not of type Distribution");
    
    Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),
                        "[disFactory] Setting up compartment distribution ", static_cast<int>(i));
    
    if (dist_i.containsElementNamed("type")) {
      Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),
                          "[disFactory] creating distribution object "
                          + Rcpp::as<std::string>(dist_i["type"]));
    }
    
    // call DistributionFactory to create Distribution variable
    auto dist_ptr = make_distribution_from_list(dist_i, component_, parameters);
    
    // call ParameterModelFactory to create ParameterModel variable
    // identify not-shape parameter, pull the parameter, and return with name
    auto modeled = extract_modeled_param(dist_i);
    Rcpp::List modeled_param = modeled.first;
    std::string modeled_param_name = modeled.second;
    
    Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),
                        "[disFactory] creating model " + modeled_param_name);
    
    Rcpp::List dp = dist_i["params"];
    
    auto model_ptr = make_param_model_from_list(component_,
                                                dp[modeled_param_name], 
                                                modeled_param_name,
                                                parameters);
    
    cds.emplace_back(std::move(dist_ptr), std::move(model_ptr));
  }
  return cds;
}
#endif
