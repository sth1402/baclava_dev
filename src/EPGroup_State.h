#ifndef EPGROUPSTATE_H
#define EPGROUPSTATE_H

#include <Rcpp.h>
#include "Pkg_types.h"

class Model;

struct EPGroup_State {

  Rcpp::List data;
  
  Rcpp::NumericVector tau;
  Rcpp::IntegerVector Z;
  
  int iP, ipsi, n, n_intvl;
  Pkg::EP endpoint_type;
  Rcpp::NumericVector endpoint_times;
  Rcpp::NumericVector age_entry;
  
  std::vector<Rcpp::NumericVector> ep_values;
  std::vector<std::vector<double>> age_screen;
  std::vector<std::vector<int>> screen_types;
  
  bool marginalize = false;
  
  // cached log_likelihood components
  mutable Rcpp::NumericVector logH, logP, logI, logB, logCP;
                      
  // constructor
  EPGroup_State(Rcpp::List data_, 
                Rcpp::NumericVector tau_, 
                Rcpp::IntegerVector Z_) : data(data_), tau(tau_), Z(Z_) {
    this->iP = Rcpp::as<int>(data["iP"]);
    this->ipsi = Rcpp::as<int>(data["ipsi"]);
    this->n = Rcpp::as<int>(data["n"]);
    
    std::string ep_type = Rcpp::as<std::string>(data["endpoint_type"]);
    if (ep_type == "clinical") this->endpoint_type = Pkg::EP::CLINICAL;
    if (ep_type == "censored") this->endpoint_type = Pkg::EP::CENSORED;
    if (ep_type == "preclinical") this->endpoint_type = Pkg::EP::PRECLINICAL;
    
    this->endpoint_times = Rcpp::as<Rcpp::NumericVector>(data["endpoint_time"]);
    this->age_entry = Rcpp::as<Rcpp::NumericVector>(data["age_entry"]);
    
    Rcpp::List endpoints_list = Rcpp::as<Rcpp::List>(data["endpoints"]);
    Rcpp::NumericVector tep = Rcpp::as<Rcpp::NumericVector>(endpoints_list["values"]);
    Rcpp::IntegerVector ep_starts = Rcpp::as<Rcpp::IntegerVector>(endpoints_list["starts"]);
    Rcpp::IntegerVector ep_ends = Rcpp::as<Rcpp::IntegerVector>(endpoints_list["ends"]);
    Rcpp::IntegerVector ep_lengths = Rcpp::as<Rcpp::IntegerVector>(endpoints_list["lengths"]);
    this->n_intvl = sum(ep_lengths) - n;
    
    for (int i = 0; i < this->n; ++i) {
      std::vector<double> tmp;
      for (int j = ep_starts[i]; j <= ep_ends[i]; ++j) {
        tmp.push_back(tep[j]);
      }
       this->ep_values.push_back(Rcpp::wrap(tmp));
    }
    
    Rcpp::List ages_screen_list = Rcpp::as<Rcpp::List>(data["ages_screen"]);
    
    if (ages_screen_list.size() == 0) {
      for (int i = 0; i < this->n; ++i) {
        std::vector<double> tmp(0);
        this->age_screen.push_back(tmp);
        std::vector<int> itmp(0);
        this->screen_types.push_back(itmp);
      }
    } else {
      Rcpp::NumericVector tap = Rcpp::as<Rcpp::NumericVector>(ages_screen_list["values"]);
      Rcpp::IntegerVector sa_starts = Rcpp::as<Rcpp::IntegerVector>(ages_screen_list["starts"]);
      Rcpp::IntegerVector sa_ends = Rcpp::as<Rcpp::IntegerVector>(ages_screen_list["ends"]);
      Rcpp::IntegerVector tst = Rcpp::as<Rcpp::IntegerVector>(ages_screen_list["types"]);
      
      for (int i = 0; i < this->n; ++i) {
        std::vector<double> tmp;
        std::vector<int> itmp;
        for (int j = sa_starts[i]; j <= sa_ends[i]; ++j) {
          tmp.push_back(tap[j]);
          itmp.push_back(tst[j]);
        }
        this->age_screen.push_back(tmp);
        this->screen_types.push_back(itmp);
      }
    }
  }
    
  void initialize_cache(const Model* model);
  
  void initialize_cache(const Model* model, const Rcpp::NumericVector& tau_, 
                        const Rcpp::IntegerVector& Z_);
  
  double calc_psi(const double tau_eval, const Model* model) const;
  Rcpp::NumericVector calc_psi(const Rcpp::NumericVector& tau_eval, const Model* model) const;
  
  Rcpp::NumericVector pfunc_ab_optimized(const std::vector<double>& dist_scales, 
                                         Pkg::Dist dist_type, 
                                         const double shape_val, const double t0, 
                                         const Rcpp::NumericVector& wgt) const;
    
  // likelihood H component
  void H_term(Rcpp::NumericVector& result,
              const Rcpp::NumericVector& tau_eval,
              const Model* model, uint8_t impact, bool fallback_to_cache = true) const;
  
  Rcpp::NumericVector H_term(const Rcpp::NumericVector& tau_eval,
                             const Model* model, uint8_t impact) const {
    Rcpp::NumericVector result(this->n, 0.0);
    H_term(result, tau_eval, model, impact, /*fallback_to_cache*/true);
    return result;
  }
  
  // likelihood P component
  void P_term(Rcpp::NumericVector& result,
              const Rcpp::NumericVector& tau_eval,
              const Rcpp::IntegerVector& indolent,
              const Model* model, uint8_t impact, bool fallback_to_cache = true) const;
  
  Rcpp::NumericVector P_term(const Rcpp::NumericVector& tau_eval,
                             const Rcpp::IntegerVector& indolent,
                             const Model* model, uint8_t impact) const {
    Rcpp::NumericVector result(this->n, 0.0);
    P_term(result, tau_eval, indolent, model, impact, /*fallback_to_cache*/true);
    return result;
  }
  
  // likelihood psi/Z component
  void I_term(Rcpp::NumericVector& result,
              const Rcpp::NumericVector& tau_eval,
              const Rcpp::IntegerVector& indolent, 
              const Model* model, uint8_t impact, bool fallback_to_cache = true) const;
  
  Rcpp::NumericVector I_term(const Rcpp::NumericVector& tau_eval,
                             const Rcpp::IntegerVector& indolent, 
                             const Model* model, uint8_t impact) const {
    Rcpp::NumericVector result(this->n);
    I_term(result, tau_eval, indolent, model, impact, /*fallback_to_cache*/true);
    return result;
  }
  
  // likelihood beta component
  void beta_term(Rcpp::NumericVector& result,
                 const Rcpp::NumericVector& tau_eval,
                 const Model* model, uint8_t impact, bool fallback_to_cache = true) const;
  
  Rcpp::NumericVector beta_term(const Rcpp::NumericVector& tau_eval,
                                const Model* model, uint8_t impact) const {
    Rcpp::NumericVector result(this->n);
    beta_term(result, tau_eval, model, impact, /*fallback_to_cache*/true);
    return result;
  }
  
  // likelihood normalization component
  // vector of length n
  void CP_term(Rcpp::NumericVector& result, const Model* model, 
               uint8_t impact, bool fallback_to_cache = true) const;
  
  Rcpp::NumericVector CP_term(const Model* model, uint8_t impact) const {
    Rcpp::NumericVector result(this->n);
    CP_term(result, model, impact, /*fallback_to_cache*/true);
    return result;
  }
  
  void update(const Model* model, uint8_t impact) {
    this->H_term(this->logH, this->tau, model, impact, /*fallback_to_cache*/false);
    this->P_term(this->logP, this->tau, this->Z, model, impact, /*fallback_to_cache*/false);
    this->I_term(this->logI, this->tau, this->Z, model, impact, /*fallback_to_cache*/false);
    this->beta_term(this->logB, this->tau, model, impact, /*fallback_to_cache*/false);
    this->CP_term(this->logCP, model, impact, /*fallback_to_cache*/false);
  }
  
  void update(const Model* model, const Rcpp::IntegerVector& newZ, uint8_t impact) {
    this->Z = newZ;
    update(model, impact);
  }
  
  void update(const Model* model, const Rcpp::NumericVector& new_tau, uint8_t impact) {
    this->tau = new_tau;
    update(model, impact);
  }
  
  void update(const Model* model, 
              const Rcpp::NumericVector& new_tau, 
              const Rcpp::IntegerVector& newZ, uint8_t impact) {
    this->tau = new_tau;
    this->Z = newZ;
    update(model, impact);
  }
  
  // partial likelihood from cached values
  Rcpp::NumericVector log_likelihood() { 
    return this->logCP + this->logP + this->logI + this->logB + this->logH;
  }

  // partial likelihood from new model; do not cache
  Rcpp::NumericVector log_likelihood(const Model* model, uint8_t impact) { 
    return this->CP_term(model, impact) +
      this->H_term(this->tau, model, impact) +
      this->P_term(this->tau, this->Z, model, impact) +
      this->I_term(this->tau, this->Z, model, impact) +
      this->beta_term(this->tau, model, impact); }

 // full likelihood from new model and latent indolence; do not cache
  Rcpp::NumericVector log_likelihood(const Model* model, 
                                     const Rcpp::IntegerVector& newZ,
                                     uint8_t impact) const { 

    return this->H_term(this->tau, model, impact) + 
      this->P_term(this->tau, newZ, model, impact) + 
      this->I_term(this->tau, newZ, model, impact) + 
      this->CP_term(model, impact) + 
      this->beta_term(this->tau, model, impact); }
  
  // full likelihood from new model and latent tau_HP; do not cache
  Rcpp::NumericVector log_likelihood(const Model* model, 
                                     const Rcpp::NumericVector& new_tau,
                                     uint8_t impact) const { 
    
    return this->H_term(new_tau, model, impact) + 
      this->P_term(new_tau, this->Z, model, impact) + 
      this->I_term(new_tau, this->Z, model, impact) + 
      this->CP_term(model, impact) + 
      this->beta_term(new_tau, model, impact); }
  
  // full likelihood from new model and latent indolence and tau_HP; do not cache
  Rcpp::NumericVector log_likelihood(const Model* model, 
                                     const Rcpp::NumericVector& new_tau,
                                     const Rcpp::IntegerVector& newZ,
                                     uint8_t impact) const { 
    
    return this->H_term(new_tau, model, impact) + 
      this->P_term(new_tau, newZ, model, impact) + 
      this->I_term(new_tau, newZ, model, impact) + 
      this->CP_term(model, impact) + 
      this->beta_term(new_tau, model, impact); }
};
#endif
