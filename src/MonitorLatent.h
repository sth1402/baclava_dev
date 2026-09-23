#ifndef MONITORLATENT_H
#define MONITORLATENT_H

#include <Rcpp.h>
#include <string>
#include <vector>

struct EPGroup_State;

class MonitorLatent {
public:
  Rcpp::NumericMatrix Z_tau;
  Rcpp::IntegerMatrix Z_indolent;
  Rcpp::LogicalMatrix tau_accept;
  Rcpp::NumericVector accept_prob;
  Rcpp::IntegerVector n_prop;
  Rcpp::IntegerVector n_acc;
  
  Rcpp::NumericVector latest_tau;
  Rcpp::IntegerVector latest_ind;
  
  int save_latent = 0;
  int save_Z = 0;
  int n_cases = 0;
  
  // Chunking / streaming
  std::string store_dir;
  int dump_every = 0; 
  int chunk_filled = 0;
  int chunk_index = 0;
  std::string file_prefix;
  
  Rcpp::IntegerVector order;
  
  // Reusable per-update buffers (avoid realloc each iteration)
  Rcpp::NumericVector tau_buf;
  Rcpp::IntegerVector indolent_buf;
  Rcpp::LogicalVector accept_buf;
  
  
  MonitorLatent() {};
  MonitorLatent(const int n_cases_, const int save_latent_,
                const bool I_marginalize,
                const std::string& store_dir_,
                const int dump_every_,
                const std::string& file_prefix_,
                const Rcpp::IntegerVector order_);
  
  
  void update(const std::vector<EPGroup_State>& groupstates,
              const Rcpp::List& accept_tau,
              const int ikept, const bool keep);
  
  void finalize();
  void flush_chunk_();
  
  
};

#endif
