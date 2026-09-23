#include <Rcpp.h>

#include "tauHP_proposal.h"
#include "ind_proposal.h"
#include "Model.h"
#include "Updates.h"
#include "Step_Factory.h"
#include "Step_MH.h"
#include "Step_Gibbs.h"
#include "Monitor.h"
#include "MonitorLatent.h"
#include "NegativeRateError.h"
#include "Verbosity.h"
#include "EPGroup_State.h"
#include "Pkg_types.h"

using namespace Rcpp;

/* Full Markov Chain Monte Carlo Procedure
 * @param data_objects List The data broken down according to modelling subsets
 *   and the endpoint type.
 * @param indolents List The current estimates for indolent aligned with data_objects.
 * @param taus List The current estimates for the age at time of 
 *   healthy -> preclinical transition aligned with data_objects.
 * @param H List the healthy compartment model specifications. 
 * @param P List the pre-clinical compartment model specifications. 
 * @param psi List the psi model specifications. 
 * @param beta List the screening sensitivity model specifications. 
 * @param t0 double The initial time.
 * @param M A positive scalar integer specifying the number of MC samples.
 * @param thin A positive scalar integer specifying the number of samples to
 *   be skipped between accepted samples.
 * @param burnin A positive scalar integer specifying the number of samples to
 *   exclude from returned posterior.
 * @param adaptive List setting for the adaptive updating of epsilons
 * @param verbose integer if 1, messages are generated showing progress
 * @param save_latent integer if 1 latent variable are kept at all kept
 *   iterations
 * @param debug integer if 1 verbose messages are generated to facilitate
 *   debugging
 * @param psi_update_type string indicating update method
 *   
 * @returns List 
 */
// [[Rcpp::export]]
List MCMC_cpp_internal(List data_objects, 
                       List indolents, 
                       List taus, 
                       List H, List P, List psi, List beta,
                       List strategies,
                       double t0, int M, int thin, int burnin,
                       int verbose, int save_latent, int debug,
                       // int use_bug, bool use_orig_prop, 
                       // bool I_marginalize_CEN, bool I_marginalize_PRE, bool I_marginalize_CLI,
                       std::string latent_dir, int latent_dump, std::string file_prefix,
                       IntegerVector latent_order) {
  
  switch (debug) {
    case 0: Verbosity::instance().setMask(NONE); break;
    case 1: Verbosity::instance().setMask(ERROR); break;
    case 2: Verbosity::instance().setMask(ERROR | INFO); break;
    case 3: Verbosity::instance().setMask(ERROR | INFO | TRACE); break;
    case 4: Verbosity::instance().setMask(ALL); break;
    default: Verbosity::instance().setMask(NONE); break;
  }

  Rcpp::RNGScope scope; // ties to R's RNG state for this call
  
  bool keep_itr;
  
  // identifying sizes for stored/returned variables
  const int N = std::max(0, M - burnin);
  const int M_thin = (N == 0) ? 0 : ( (N + thin - 1) / thin );
  if (M_thin < 1) stop("M_thin is 0; verify inputs");
  
  int ikept = -1;
  
  // extract the number of participants in each subset of the data
  int total_cases = 0;
  for(R_xlen_t k = 0; k < data_objects.length(); ++k) {
    List dobj = data_objects[k];
    total_cases += as<int>(dobj["n"]);
  }
  Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),  
                "Total # participants ", total_cases);

  // initialize latent monitoring object
  MonitorLatent monitor_latent(total_cases, save_latent, 
                               false,
                               latent_dir, latent_dump, file_prefix, latent_order);
  
  Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),  
                "preparing model");

  // prepare model
  Model model = Model();
  try {
    model = Model(H, P, psi, beta, t0);
  } catch (const NegativeRateError& error) {
    stop("could not create model");
  }
  model.I_marginalize_CEN = false;
  model.I_marginalize_PRE = false;
  model.I_marginalize_CLI = false;

  std::vector<Step_MH> MH_steps;
  std::vector<Step_Gibbs> Gibbs_steps;
  
  try {
    init_Step_MHs(strategies, model.parameters, MH_steps);
  } catch (const std::exception& e) {
    stop("failed to create MH strategies: %s", e.what());
  } catch (...) {
    stop("failed to create MH strategies with unknown exception");
  }
  
  try {
    init_Step_Gibbs(strategies, model.parameters, Gibbs_steps);
  } catch (const std::exception& e) {
    stop("failed to create Gibbs strategies: %s", e.what());
  } catch (...) {
    stop("failed to create Gibbs strategies with unknown exception");
  }
  
  struct UpdateStep {
    Pkg::Step type;
    int group_index; 
  };
  
  std::vector<UpdateStep> procedure;
  int MH_idx = 0;
  int GB_idx = 0;
  List blocks = strategies["blocks"];
  int nb = strategies["n_blocks"];
  for (int i = 0; i < nb; ++i) {
    List block = as<List>(blocks[i]);
    std::string type = as<std::string>(block["alt_meta_for"]);
    if (type == "MH") {
      UpdateStep us;
      us.type = Pkg::Step::MH;
      us.group_index = MH_idx;
      procedure.push_back(us);
      MH_idx++;
    } else if (type == "MH_with_Z") {
      UpdateStep us;
      us.type = Pkg::Step::MH_with_Z;
      us.group_index = MH_idx;
      procedure.push_back(us);
      MH_idx++;
    } else if (type == "Gibbs") {
      UpdateStep us;
      us.type = Pkg::Step::Gibbs;
      us.group_index = GB_idx;
      procedure.push_back(us);
      GB_idx++;
    }
  }

  Verbosity::instance().log(TRACE, "model created ");

  // initialize monitors
  std::vector<Monitor> MH_monitors;
  for (auto& mh : MH_steps) MH_monitors.push_back(Monitor(mh, M_thin));
  std::vector<Monitor> Gibbs_monitors;
  for (auto& gibbs : Gibbs_steps) Gibbs_monitors.push_back(Monitor(gibbs, M_thin));
  
  // initialize tracing of acceptance of tau at the individual level
  List accept_tau(data_objects.size());

  // groupstates are an internal convenience function to cache likelihood
  // calculations and to store data information conveniently
  std::vector<EPGroup_State> groupstates;
  
  if (taus.size() != data_objects.size()) {
    // if not a continuation step, must initialize latent variables tau and Z
    // With the tau dependent rate models, need to provide initial 
    // guesses for tau
    NumericVector prob_tau;
    NumericVector prob_indolents;
    taus = List(data_objects.size());
    indolents = List(data_objects.size());
    
    for (int i = 0; i < data_objects.size(); i++) {

      List dobj = as<List>(data_objects[i]);
      
      std::string endpoint_type = as<std::string>(dobj["endpoint_type"]);
      NumericVector endpoint_time = as<NumericVector>(dobj["endpoint_time"]);
      
      int n = as<int>(dobj["n"]);

      if (endpoint_type == "censored") {
        NumericVector tmp = endpoint_time + runif(n, -1.0, 1.0);
        tmp = ifelse(tmp < model.t0, model.t0 + 0.1, tmp);
        taus[i] = tmp;
        
        IntegerVector indolent(n, 0);
        indolents[i] = indolent;

      } else if (endpoint_type == "preclinical") {
        NumericVector tmp = endpoint_time - runif(n, 0.5, 1.0);
        tmp = ifelse(tmp < model.t0, model.t0 + 0.1, tmp);
        taus[i] = tmp;
        
        IntegerVector indolent(n, 0);
        indolents[i] = indolent;

      } else {
        // for all clinical, sample tau, but indolent is 0
        NumericVector tmp = endpoint_time - runif(n, 0.5, 5.0);
        taus[i] = ifelse(tmp < model.t0, model.t0 + 0.1, tmp);
        indolents[i] = IntegerVector(n, 0);
      }
      EPGroup_State gs(data_objects[i], taus[i], indolents[i]);
      gs.initialize_cache(&model);
      groupstates.push_back(gs);
    }
    
    Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),  
                        "cycling through seeding of tau and Z");
    // cycle through a few times to stabilize initial values
    for (int i = 0; i < 10; i++) {
      for (int j = 0; j < groupstates.size(); ++j) {
        auto& gs = groupstates[j];
        NumericVector wgt = compute_prob_tau_prep(gs, &model);
        prob_tau = compute_prob_tau(gs, gs.tau, wgt, &model);
        gs.tau = rprop_tau(gs, gs.tau, prob_tau, &model);
        gs.update(&model, static_cast<uint8_t>(Pkg::Components::TAU));
        
        prob_indolents = compute_prob_indolent(gs, gs.tau, &model);
        gs.Z = rprop_indolent(gs.endpoint_type, prob_indolents);
        gs.update(&model, static_cast<uint8_t>(Pkg::Components::Z));
      }
    }
    
    if (debug > 0) {
      for (int i = 0; i < groupstates.size(); i++) {
        NumericVector tmp2 = groupstates[i].tau;
        IntegerVector itmp = groupstates[i].Z;
        int n = groupstates[i].n;
        Verbosity::instance().log(TRACE,  
                            "after compute_prob_indolent_List " + std::to_string(i) + " " 
                                    + std::to_string(mean(tmp2)) + " " 
                                     + std::to_string(mean(itmp)) + " " 
                                    + std::to_string(n));
        int cnt = 0;
        for (int j = 0; j < groupstates[i].n; ++j)
          if (groupstates[i].tau[j] > groupstates[i].endpoint_times[j]) cnt++;
        Rcout << cnt << std::endl;
      }
    }
  } else {
    // if a continuation step, initialize group states based on provided latents
    for (int i = 0; i < data_objects.size(); ++i) {
      EPGroup_State gs(data_objects[i], taus[i], indolents[i]);
      gs.initialize_cache(&model);
      groupstates.push_back(gs);
    }
  }
  
  // progression screen prints
  int m_mod = M;
  if (verbose == 1) {
    m_mod = floor(M * 0.02);
    if (m_mod > 0) {
      Rcout << "0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%" << std::endl;
      Rcout << "[----|----|----|----|----|----|----|----|----|----|" << std::endl;
    }
  }
  
  for (int m = 0; m < M; ++m) {

    Rcpp::checkUserInterrupt();
    if (verbose == 1 && m_mod > 0 && m % m_mod == 0) Rcout << "*" ;
    
    Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),  
                  "Iteration", m);
    
    // track if iteration will be kept
    keep_itr = m >= burnin && ((m - burnin) % thin == 0);
    if (keep_itr) ikept += 1;
    
    for (auto step : procedure) {
      if (step.type == Pkg::Step::MH) {
        try {
          Updates::MH(groupstates, &model, MH_steps[step.group_index], m);
          if (m >= burnin) 
            MH_monitors[step.group_index].update(MH_steps[step.group_index], ikept, keep_itr);
        } catch (const std::exception& e) {
          stop("failed to update MH: %s", e.what());
        } catch (...) {
          stop("failed to update MH with unknown exception");
        }
      } else if (step.type == Pkg::Step::Gibbs) {
        try {
          Updates::beta(groupstates, &model, Gibbs_steps[step.group_index]);
          if (m >= burnin) 
            Gibbs_monitors[step.group_index].update(Gibbs_steps[step.group_index], ikept, keep_itr);
        } catch (const std::exception& e) {
          stop("failed to update Gibbs beta: %s", e.what());
        } catch (...) {
          stop("failed to update Gibbs beta with unknown exception");
        }
      }
    }
    
    // update tau and monitor accepted latent variables Z and tau
    try {
      Updates::tau(groupstates, accept_tau, &model);
      monitor_latent.update(groupstates, accept_tau, ikept, keep_itr);
    } catch (const std::exception& e) {
      stop("failed to update tau: %s", e.what());
    } catch (...) {
      stop("failed to update tau with unknown exception");
    }
    
    for (auto step : procedure) {
      if (step.type == Pkg::Step::MH_with_Z) {
        try {
          Updates::psiZ(groupstates, &model, MH_steps[step.group_index], m);
          if (m >= burnin) 
            MH_monitors[step.group_index].update(MH_steps[step.group_index], ikept, keep_itr);
        } catch (const std::exception& e) {
          stop("failed to update psi: %s", e.what());
        } catch (...) {
          stop("failed to update psi with unknown exception");
        }
      }
    }
    
  } // end runtime

  if (verbose == 1) {
    Rcout << std::endl << std::endl;
    Rcout << "sampler completed" << std::endl;
  }

  // extracting evolution of update variable s
  List adaptive_s;
  try {
    Verbosity::instance().log(TRACE,  
                  "constructing adaptive history");
    for (auto& mh : MH_steps) {
      CharacterVector names = mh.param_names;
      std::string nm = "";
      for (auto st : names) nm = nm + "/" + st;
      if (!nm.empty()) nm.erase(0, 1);
      adaptive_s[nm] = mh.getEvolution();
    }
  } catch (const std::exception& e) {
    stop("failed to extract evolution of s: %s", e.what());
  } catch (...) {
    stop("failed to extract evolution of s with unknown exception");
  }
  
  // extracting update settings (s, L, d)
  List k_epsilon;
  try {
    Verbosity::instance().log(TRACE, "extracting update settings");
    
    for (auto& mh : MH_steps) {
      CharacterVector names = mh.param_names;
      std::string nm = "";
      for (auto st : names) nm = nm + "/" + st;
      if (!nm.empty()) nm.erase(0, 1);
      k_epsilon[nm] = mh.extract_epsilon();
    }
  } catch (const std::exception& e) {
    stop("failed to extract update settings: %s", e.what());
  } catch (...) {
    stop("failed to extract update settings with unknown exception");
  }
  
  // output
  Verbosity::instance().log(TRACE, "creating returned list");
  
  // kept posteriors
  List kTheta;
  List kACCEPT;
  List acceptance_prob;
  for (auto& mh : MH_monitors) {
    CharacterVector names = mh.param_names;
    std::string nm = "";
    for (auto st : names) nm = nm + "/" + st;
    if (!nm.empty()) nm.erase(0, 1);
    kTheta[nm] = mh.posterior;
    kACCEPT[nm] = mh.kACCEPT;
    mh.update_acceptance_prob();
    acceptance_prob[nm] = mh.acceptance_prob;
  }
  
  for (auto& gb : Gibbs_monitors) {
    CharacterVector names = gb.param_names;
    std::string nm = "";
    for (auto st : names) nm = nm + "/" + st;
    if (!nm.empty()) nm.erase(0, 1);
    kTheta[nm] = gb.posterior;
  }
  
  List latent;
  if (save_latent == 1) {
    monitor_latent.finalize();
    latent["latent_dir"] = monitor_latent.store_dir;
    latent["file_prefix"] = monitor_latent.file_prefix;
    latent["latent_dump"] = monitor_latent.dump_every;
    latent["tau_accept_rate"] = monitor_latent.accept_prob[monitor_latent.order];
    latent["last_iter_tau"] = monitor_latent.latest_tau[monitor_latent.order];
    latent["last_iter_indolent"] = monitor_latent.latest_ind[monitor_latent.order];
  } else {
    NumericVector tmp = monitor_latent.Z_tau(0, _);
    latent["tau"] = tmp[monitor_latent.order];
    tmp = monitor_latent.tau_accept(0, _);
    latent["tau_accept"] = tmp[monitor_latent.order];
    latent["tau_accept_rate"] = monitor_latent.accept_prob[monitor_latent.order];
    if (monitor_latent.save_Z) {
      IntegerVector itmp = monitor_latent.Z_indolent(0, _);
      latent["indolent"] = itmp[monitor_latent.order];
    }
  }
  
  return List::create(Named("theta") = kTheta, 
                      Named("latent") = latent,
                      Named("accept") = kACCEPT,
                      Named("acceptance_prob") = acceptance_prob,
                      Named("epsilon") = k_epsilon,
                      Named("adaptive") = adaptive_s,
                      Named("ikept") = ikept + 1); 
}
