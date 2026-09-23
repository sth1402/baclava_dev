#ifndef MODEL_H
#define MODEL_H

#include <Rcpp.h>
#include "Compartment_Factory.h"

class Compartment;
class StateParameter;

class Model {
public:
  
  std::vector<std::shared_ptr<StateParameter>> parameters;
  
  std::vector<Compartment> H;
  std::vector<Compartment> P;
  std::vector<Compartment> psi;
  std::vector<Compartment> beta;
  
  double t0;
  
  Model() {};
  Model(Rcpp::List H_, Rcpp::List P_, Rcpp::List psi_, Rcpp::List beta_, 
        double t0);
  
  bool I_marginalize_CEN = false;
  bool I_marginalize_PRE = false;
  bool I_marginalize_CLI = false;
  
private:
  
  inline std::vector<Compartment> init_pack(
      const Rcpp::List& comp, Pkg::Components component_,
      std::vector<std::shared_ptr<StateParameter>>& parameters) {
    
    // ensure that comp is appropriately formatted
    if (!comp.containsElementNamed("meta_for"))
      Rcpp::stop("Model constructor: provided compartment list invalid");
    
    if (Rcpp::as<std::string>(comp["meta_for"]) != "Compartment")
      Rcpp::stop("init_pack provided list is not of type Compartment");
    
    if (!comp.containsElementNamed("distributions"))
      Rcpp::stop("comp is not properly formatted");
    
    // build Compartment component
    
    return build_compartment(comp["distributions"], component_, parameters);
    
  }
  
};

#endif
