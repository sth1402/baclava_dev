#include <Rcpp.h>
#include "Step_Factory.h"
#include "Compartment.h"
#include "Step_MH.h"
#include "Step_Gibbs.h"
#include "StateParameter.h"

using namespace Rcpp;

ParamIndex index_parameters_by_name(const std::vector<Compartment>& cds) { 
  ParamIndex idx; 
  
  for (int c = 0; c < static_cast<int>(cds.size()); ++c) { 
    const auto& cd = cds[c]; 
    
    // distribution shape, if present 
    if (cd.distribution && cd.distribution->shape) { 
      const auto& sp = cd.distribution->shape; 
      if (sp) idx[sp->tag].emplace_back(sp, c);
    } 
    
    // modeled parameter(s) 
    if (cd.parameterModel) { 
      for (const auto& p : cd.parameterModel->parameters) { 
        if (!p) continue; 
        idx[p->tag].emplace_back(p, c); 
      } 
    } 
  } 
  return idx; 
}

/*
 * BlockStrategies::toList = function() {
 *   list(
 *     meta_for = "BlockStrategies",
 *     n_blocks = length(self$blocks),
 *     blocks = lapply(self$blocks, function(b) b$toList())
 *   )
 * }
 */
void setup_compartment_updates_MH(
    const List block_strategies,
    const std::vector<std::shared_ptr<StateParameter>>& parameters,
    std::vector<Step_MH>& MH_steps) {
  
  // verify that the block_strategies is appropriately defined
  
  if (!block_strategies.containsElementNamed("meta_for"))
    stop("setup_compartment_updates_MH: provided parameter list invalid");
  
  if (as<std::string>(block_strategies["meta_for"]) != "Strategies")
    stop("setup_compartment_updates_MH: provided list is not of type Strategies");
  
  if (!block_strategies.containsElementNamed("n_blocks") ||
      !block_strategies.containsElementNamed("blocks"))
      stop("block_strategies is not properly formatted");
  
  // number of block strategies for this compartment
  const int n_blocks = as<int>(block_strategies["n_blocks"]);
  
  if (n_blocks == 0) return;

  List blocks = block_strategies["blocks"];
  
  // for each block strategy, create an Step_MH object
  for (int b = 0; b < n_blocks; ++b) {
    
    List block = as<List>(blocks[b]);
    
    // ensure that block object is appropriately defined
    
    if (!block.containsElementNamed("meta_for"))
      stop("setup_compartment_updates: provided Block list invalid");
    
    if (as<std::string>(block["meta_for"]) != "MH" && as<std::string>(block["meta_for"]) != "Gibbs")
      stop("setup_compartment_updates: provided list is not of type MH or Gibbs");
    
    if (as<std::string>(block["meta_for"]) != "MH") continue;
    
    if (!block.containsElementNamed("params"))
        stop("block is not properly formatted");
    
    // extract the StateParameter tags included in the MH block
    CharacterVector param_names = block["params"];
    
    if (param_names.size() == 0) stop("cannot provide an empty parameter list in block");
    
    // new vector of pointers to StateParameters specific to this update block
    std::vector<std::shared_ptr<StateParameter>> group_params;
    group_params.reserve(param_names.size());
    
    for (R_xlen_t i = 0; i < param_names.size(); ++i) {
      const std::string nm = as<std::string>(param_names[i]);
      
      bool already_there = false;
      for (auto& p : parameters) {
        if (p->tag == nm) {
          already_there = true;
          group_params.push_back(p);
          break;
        }
      }
      
      if (!already_there) stop("strategy parameter not found in model" + nm);
    }    
    if (group_params.empty()) stop("Block matched no updatable parameters");
    
    // Create the MH group
    MH_steps.emplace_back(group_params, param_names);
    // initialize corr_block
    MH_steps.back().init_corr_block(block);
  }
  
}


void setup_compartment_updates_Gibbs(
    const List block_strategies,
    const std::vector<std::shared_ptr<StateParameter>>& parameters,
    std::vector<Step_Gibbs>& Gibbs_groups) {
  
  // verify that the block_strategies is appropriately defined
  
  if (!block_strategies.containsElementNamed("meta_for"))
    stop("setup_compartment_updates_MH: provided parameter list invalid");
  
  if (as<std::string>(block_strategies["meta_for"]) != "Strategies")
    stop("setup_compartment_updates_MH: provided list is not of type Strategies");
  
  if (!block_strategies.containsElementNamed("n_blocks") ||
      !block_strategies.containsElementNamed("blocks"))
      stop("block_strategies is not properly formatted");
  
  // number of block strategies for this compartment
  const int n_blocks = as<int>(block_strategies["n_blocks"]);
  
  if (n_blocks == 0) return;
  
  List blocks = block_strategies["blocks"];
  
  // for each block strategy, create an Step_MH object
  for (int b = 0; b < n_blocks; ++b) {
    
    List block = as<List>(blocks[b]);
    
    // ensure that block object is appropriately defined
    
    if (!block.containsElementNamed("meta_for"))
      stop("setup_compartment_updates: provided Block list invalid");
    
    if (as<std::string>(block["meta_for"]) != "MH" && as<std::string>(block["meta_for"]) != "Gibbs")
      stop("setup_compartment_updates: provided list is not of type MH or Gibbs");
    
    if (as<std::string>(block["meta_for"]) != "Gibbs") continue;
    
    if (!block.containsElementNamed("params"))
      stop("block is not properly formatted");
    
    // extract the StateParameter tags included in the MH block
    CharacterVector param_names = block["params"];
    
    if (param_names.size() == 0) stop("cannot provide an empty parameter list in block");
    
    // new vector of pointers to StateParameters specific to this update block
    std::vector<std::shared_ptr<StateParameter>> group_params;
    group_params.reserve(param_names.size());
    
    for (R_xlen_t i = 0; i < param_names.size(); ++i) {
      const std::string nm = as<std::string>(param_names[i]);
      
      bool already_there = false;
      for (auto& p : parameters) {
        if (p->tag == nm) {
          already_there = true;
          group_params.push_back(p);
          break;
        }
      }
      
      if (!already_there) stop("strategy parameter not found in model");
    }    
    if (group_params.empty()) stop("Block matched no updatable parameters");
    
    // Create the MH group
    Gibbs_groups.emplace_back(group_params, param_names);
  }
}