#ifndef CORRBLOCK_H
#define CORRBLOCK_H

#include <Rcpp.h>
#include <vector>
#include <memory>
#include <limits>
#include "EpsAdapter.h"
#include "CovAdapter.h"
#include "Verbosity.h"

class StateParameter;
class Proposal_State;

// The proposal procedure for a single MH update step
/*
 * CorrBlock
 * 
 * Coordinates step-size and covariance adaptation for an MCMC parameter block.
 *   - Uses EpsAdapter (dual averaging) to tune global scale s.
 *   - Uses CovAdapter (Welford + Cholesky) to adapt correlation structure.
 */

struct CorrBlock {
  
  CorrBlock() {};
  
  /*
   *  Block$toList = function() {
   *    list(
   *      "meta_for" = "Block",
   *      "params" = self$params,
   *      "s0" = self$s0,
   *      "adaptive" = self$adaptive,
   *.     "L" = self$L,
   *      "D" = self$D
   *    )
   * }
   */
  explicit CorrBlock(const Rcpp::List& block_list)  {
    
    Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),
                  "[corrBlock] initializing from list");
    
    if (!block_list.containsElementNamed("meta_for"))
      Rcpp::stop("CorrBlock Constructor: provided Block list invalid");
    
    if (Rcpp::as<std::string>(block_list["meta_for"]) != "MH")
      Rcpp::stop("CorrBlock Constructor: provided list is not of type MH");
    
    if (!block_list.containsElementNamed("s0") ||
        !block_list.containsElementNamed("adaptive") ||
        !block_list.containsElementNamed("L") ||
        !block_list.containsElementNamed("D"))
        Rcpp::stop("block_list is not properly formatted");
    
    try{
      this->s = Rcpp::as<double>(block_list["s0"]);
      this->L = Rcpp::as<Rcpp::NumericMatrix>(block_list["L"]);
      this->D = Rcpp::as<Rcpp::NumericVector>(block_list["D"]);
      this->d = L.nrow();

      // create covariance and eps adapters
      const Rcpp::List adapt_list = Rcpp::as<Rcpp::List>(block_list["adaptive"]);
      
      this->cov_adapter = CovAdapter(adapt_list);
      this->eps_adapter = EpsAdapter(this->s, adapt_list, this->d);
      
      if (adapt_list.containsElementNamed("warmup")) {
        this->warmup = Rcpp::as<int>(adapt_list["warmup"]);
      } else {
        this->warmup = 0;
      }
      
      if (adapt_list.containsElementNamed("warmup_cov")) {
        this->warmup_cov = Rcpp::as<int>(adapt_list["warmup_cov"]);
        if (this->warmup_cov == 0) this->warmup_cov = -1;
      } else {
        this->warmup_cov = -1; // never update D/L
      }
      
      // store local flags for what components are to be learned
      std::string learning = Rcpp::as<std::string>(adapt_list["learn"]);
      
      if (!(learning == "none" || learning == "s" || learning == "sD" || learning == "sDL"))
        Rcpp::stop("learn must be one of {s, sD, sDL, none}");

      this->learn_L = learning == "sDL";
      this->learn_D = learning == "sDL" || learning == "sD";
      this->learn_s = learning != "none";
      this->frozen = !this->learn_s;
      
      if (this->warmup > 0) {
        if (learning == "none") Rcpp::stop("learn must be s, sD, or sDL to trigger adaptive procedure");
        this->s_evolution.reserve(this->warmup);
        if (this->learn_D) this->d_evolution.reserve(this->warmup);
        if (this->learn_L) this->l_evolution.reserve(this->warmup);
      } else {
        if (learning != "none") Rcpp::stop("warmup must be > 0 for learn = s, sD, or sDL");
        this->frozen = true;
        this->cov_adapter.freeze();
        this->eps_adapter.freeze();
      }
      
      if (this->warmup_cov < 0) this->cov_adapter.freeze();

    } catch (const std::exception& e) {
      Rcpp::stop("CorrBlock constructor: missing or malformed parameter list: %s", e.what());
    }
    
    std::string msg = "[corrBlock] adapt=" + std::to_string(!this->frozen);
    Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO), msg);

  }
  
  void adapt(const double a_hat,
              const std::vector<std::shared_ptr<StateParameter>>& parameters,
              const int m);
  
  // proposal of new values on the update scale given current values 
  // on the update (working) scale
  Rcpp::NumericVector propose(const Rcpp::NumericVector& w_old_block,
                              const Proposal_State& params) const;
  
  // return current update settings
  Rcpp::List getFinalValues() const {
    return Rcpp::List::create(
      Rcpp::Named("s") = this->s,
      Rcpp::Named("L") = this->L,
      Rcpp::Named("D") = this->D,
      Rcpp::Named("d") = this->d
    );
  }
  
  Rcpp::List getAdaptivePhase() const { 
    
    if (!this->learn_s) return Rcpp::List::create(
        Rcpp::Named("s") = NA_REAL,
        Rcpp::Named("L") = NA_REAL,
        Rcpp::Named("D") = NA_REAL,
        Rcpp::Named("target_delta") = NA_REAL
    );
    
    Rcpp::List adaptive_phase;
    adaptive_phase["s"] = this->getEvolution();
    if (this->learn_D) adaptive_phase["D"] = this->getDEvolution();
    if (this->learn_L) adaptive_phase["L"] = this->getLEvolution();
    adaptive_phase["target_delta"] = this->eps_adapter.getDelta();
    return adaptive_phase;
  }
  
  // return the values of D during the adaptive period and the
  // target acceptance rate
  CorrBlock(const CorrBlock&) = default;
  CorrBlock& operator=(const CorrBlock&) = default;
  CorrBlock(CorrBlock&&) = default;
  CorrBlock& operator=(CorrBlock&&) = default;
  
  bool isFrozen() { return this->frozen; }
  
private:
  
  double s;  // scalar step size
  Rcpp::NumericMatrix L; // correlation matrix
  Rcpp::NumericVector D; // shape vector

  bool frozen;
  int warmup;
  int warmup_cov;
  
  // number of parameters in the MH update block
  int d;
  
  // s adaptor
  EpsAdapter eps_adapter;
  CovAdapter cov_adapter;
  
  bool learn_D = true;  // true=adapt D with s
  bool learn_L = true;  // true=adapt L with s and use correlated proprosals
  bool learn_s = true;  // true=adapt
  
  // evolution logs
  std::vector<double> s_evolution; // s(t) = exp(log_eps) during warmup
  std::vector<double> sbar_evolution; // exp(log_eps_bar)
  std::vector<double> ahat_evolution; // a_hat fed to s adaptation
  std::vector<std::vector<double>> d_evolution;
  std::vector<std::vector<double>> l_evolution;
  
  void adapter_track(double a_hat);
  
  // return the values of s during the adaptive period
  Rcpp::NumericVector getEvolution() const { 
    return Rcpp::wrap(s_evolution);
  }
  
  // return the values of D during the adaptive period
  Rcpp::NumericMatrix getDEvolution() const {
    if (d_evolution.empty() || !learn_D) return Rcpp::NumericMatrix(0, 0);
    
    const int n_iter = static_cast<int>(d_evolution.size());
    const int dim_d  = static_cast<int>(d_evolution[0].size());
    
    Rcpp::NumericMatrix tmp(n_iter, dim_d);
    
    for (int i = 0; i < n_iter; ++i) {
      const std::vector<double>& d_i = d_evolution[i];
      const int len = std::min(dim_d, static_cast<int>(d_i.size()));
      for (int j = 0; j < len; ++j) tmp(i, j) = d_i[j];
    }
    
    return tmp;
  }
  
  // return the upper triangle values of L during the adaptive period
  Rcpp::NumericMatrix getLEvolution() const {
    if (l_evolution.empty() || !learn_L) return Rcpp::NumericMatrix(0, 0);
    
    const int n_iter = static_cast<int>(l_evolution.size());
    const int dim_d  = static_cast<int>(l_evolution[0].size());
    
    Rcpp::NumericMatrix tmp(n_iter, dim_d);
    
    for (int i = 0; i < n_iter; ++i) {
      const std::vector<double>& l_i = l_evolution[i];
      const int len = std::min(dim_d, static_cast<int>(l_i.size()));
      for (int j = 0; j < len; ++j) tmp(i, j) = l_i[j];
    }
    
    return tmp;
  }
  
};

#endif // CORRBLOCK_H
