#include <algorithm>
#include <Rcpp.h>
#include <string>
#include <map>
#include "Compartment_Factory.h"
#include "Model.h"
#include "Pkg_types.h"
#include "Verbosity.h"

using namespace Rcpp;

/*
 * Compartment$toList = function() {
 *   base <- list(
 *     "subset_var" = ifelse(self$composition == "stratified", self$subset_var, NA_character_),
 *     "labels" = names(self$distributions),
 *     "distributions" = lapply(self$distributions, function(x) x$toList()))
 *   names(base$distributions) <- names(self$distributions)
 *  base
 * }
 */


Model::Model(List H_, List P_, List psi_, List beta_, double t0_) {
  
  if (Rcpp::as<std::string>(P_["meta_for"]) != "Compartment")
    Rcpp::stop("Model constructor: provided list is not of type Compartment");
  
  // risk onset age
  this->t0 = t0_;

  // H compartment
  Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),  
                "H compartment in Model constructor");
  
  try {
    H = init_pack(H_, Pkg::Components::H, parameters);
  } catch (const std::exception& e) {
    stop("failed to create H: %s", e.what());
  } catch (...) {
    stop("failed to create H with unknown exception");
  }
  
  // P compartment
  Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),  
                "P compartment in Model constructor");
  
  try {
    P = init_pack(P_, Pkg::Components::P, parameters);
  } catch (const std::exception& e) {
    stop("failed to create P: %s", e.what());
  } catch (...) {
    stop("failed to create P with unknown exception");
  }
  
  // psi compartment
  Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),  
                "psi compartment in Model constructor");
  
  try {
    psi = init_pack(psi_, Pkg::Components::PSI, parameters);
  } catch (const std::exception& e) {
    stop("failed to create psi: %s", e.what());
  } catch (...) {
    stop("failed to create psi with unknown exception");
  }
  
  // beta compartment
  Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),  
                "beta compartment in Model constructor");
  
  try {
    beta = init_pack(beta_, Pkg::Components::BETA, parameters);
  } catch (const std::exception& e) {
    stop("failed to create beta: %s", e.what());
  } catch (...) {
    stop("failed to create beta with unknown exception");
  }
  
}