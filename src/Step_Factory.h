#ifndef STEP_FACTORY_H
#define STEP_FACTORY_H

#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <vector>
#include <string>

class StateParameter;
class Compartment;
class Step_MH;
class Step_Gibbs;

// name -> StateParameter ptr
using ParamIndex = std::unordered_map<std::string, std::vector<std::pair<std::shared_ptr<StateParameter>, int>>>; 

// Build an index from compartment distribution contents (shape + modeled params)
ParamIndex index_parameters_by_name(const std::vector<Compartment>& cds);

// Build MH groups based on an R list of block strategies
void setup_compartment_updates_MH(
    const Rcpp::List block_strategies,
    const std::vector<std::shared_ptr<StateParameter>>& parameters,
    std::vector<Step_MH>& MH_groups
);

void setup_compartment_updates_Gibbs(
    const Rcpp::List block_strategies,
    const std::vector<std::shared_ptr<StateParameter>>& parameters,
    std::vector<Step_Gibbs>& Gibbs_groups
);

inline void init_Step_MHs(
        const Rcpp::List& strategies,
        const std::vector<std::shared_ptr<StateParameter>>& parameters,
        std::vector<Step_MH>& MH_groups) {
    
    // ensure that comp is appropriately formatted
    if (!strategies.containsElementNamed("meta_for"))
        Rcpp::stop("Step_MH constructor: provided compartment list invalid");
    
    if (Rcpp::as<std::string>(strategies["meta_for"]) != "Strategies")
        Rcpp::stop("init_Step_MHs provided list is not of type Strategies");
    
    if (!strategies.containsElementNamed("blocks"))
        Rcpp::stop("strategies is not properly formatted");
    
    // build Step_MHs component
    setup_compartment_updates_MH(strategies, parameters, MH_groups);
}

inline void init_Step_Gibbs(
        const Rcpp::List& strategies,
        const std::vector<std::shared_ptr<StateParameter>>& parameters,
        std::vector<Step_Gibbs>& Gibbs_groups) {
    
    // ensure that comp is appropriately formatted
    if (!strategies.containsElementNamed("meta_for"))
        Rcpp::stop("Step_MH constructor: provided compartment list invalid");
    
    if (Rcpp::as<std::string>(strategies["meta_for"]) != "Strategies")
        Rcpp::stop("init_Step_MHs provided list is not of type Strategies");
    
    if (!strategies.containsElementNamed("blocks"))
        Rcpp::stop("strategies is not properly formatted");
    
    // build Step_MHs component
    setup_compartment_updates_Gibbs(strategies, parameters, Gibbs_groups);
}


#endif
