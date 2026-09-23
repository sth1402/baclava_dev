#ifndef DISTRIBUTION_FACTORY_H
#define DISTRIBUTION_FACTORY_H
#include "Distribution_Weibull.h"
#include "Distribution_Gamma.h"
#include "Distribution_Bernoulli.h"
#include "Verbosity.h"
#include "Pkg_types.h"

// Convenience functions to create Distribution objects

/*
 *  Distribution$toList = function() {
 *    params <- lapply(self$parameters, function(x) x$toList())
 *    names(params) <- names(self$parameters)
 *    list("meta_for" = "Distribution", 
 *         "type" = self$type,
 *         "params" = params)
 *  }
 */
inline std::shared_ptr<Distribution>
  make_distribution_from_list(
    const Rcpp::List& d, Pkg::Components component_,
    std::vector<std::shared_ptr<StateParameter>>& parameters) {

    if (!d.containsElementNamed("meta_for"))
      Rcpp::stop("Distribution Factory: provided distribution list invalid");
    
    if (Rcpp::as<std::string>(d["meta_for"]) != "Distribution")
      Rcpp::stop("Distribution Factory: provided list is not of type Distribution");
    
    if (!d.containsElementNamed("type"))
      Rcpp::stop("distribution does not contain element type");
    const std::string distribution_type = Rcpp::as<std::string>(d["type"]);
    
    if (distribution_type == "dweibull") {
      return std::make_shared<Distribution_Weibull>(d, component_, parameters);
    } else if (distribution_type == "dgamma") {
      return std::make_shared<Distribution_Gamma>(d, component_, parameters);
    } else if (distribution_type == "dbinom") {
      return std::make_shared<Distribution_Bernoulli>(d, parameters);
    }
    throw std::invalid_argument("Distribution type not recognized: " + distribution_type);
  }
#endif
