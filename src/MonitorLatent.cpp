#include <Rcpp.h>
#include "MonitorLatent.h"
#include "Verbosity.h"
#include "EPGroup_State.h"

using namespace Rcpp;

// [[plugins(cpp11)]]

MonitorLatent::MonitorLatent(const int n_cases_,
                             const int save_latent_,
                             const bool I_marginalize,
                             const std::string& store_dir_,
                             const int dump_every_,
                             const std::string& file_prefix_,
                             const IntegerVector order_)
  : save_latent(save_latent_),
    save_Z(!I_marginalize),
    n_cases(n_cases_),
    store_dir(store_dir_),
    dump_every(dump_every_),
    file_prefix(file_prefix_),
    order(order_) {
  
  accept_prob = NumericVector(n_cases, 0.0);
  n_prop      = IntegerVector(n_cases, 0);
  n_acc       = IntegerVector(n_cases, 0);
  
  tau_buf     = NumericVector(n_cases);
  indolent_buf= IntegerVector(n_cases);
  accept_buf  = LogicalVector(n_cases);
  
  if (save_latent == 1) {
    if (dump_every <= 0) stop("dump_every must be > 0 when save_latent=1");
    if (store_dir.empty()) stop("store_dir must be non-empty when save_latent=1");
    
    Z_tau = NumericMatrix(dump_every, n_cases);
    tau_accept = LogicalMatrix(dump_every, n_cases);
    if (save_Z) Z_indolent = IntegerMatrix(dump_every, n_cases);
    
    chunk_filled = 0;
    chunk_index  = 0;
    
  } else {
    Z_tau      = NumericMatrix(1, n_cases);
    tau_accept = LogicalMatrix(1, n_cases);
    if (save_Z) Z_indolent = IntegerMatrix(1, n_cases);
  }
  
  latest_tau = NumericVector(n_cases);
  latest_ind = IntegerVector(n_cases);
}

void MonitorLatent::flush_chunk_() {
  if (save_latent != 1) return;
  if (chunk_filled == 0) return;
  
  IntegerVector rows = seq(0, chunk_filled - 1);
  
  Rcpp::NumericMatrix Z_tau_out(chunk_filled, n_cases);
  Rcpp::LogicalMatrix tau_accept_out(chunk_filled, n_cases);
  Rcpp::IntegerMatrix Z_indolent_out;
  
  if (save_Z) {
    Z_indolent_out = Rcpp::IntegerMatrix(chunk_filled, n_cases);
  }
  
  for (int i = 0; i < chunk_filled; ++i) {
    for (int j = 0; j < n_cases; ++j) {
      Z_tau_out(i, j) = Z_tau(i, order[j]);
      tau_accept_out(i, j) = tau_accept(i, order[j]);
      if (save_Z) Z_indolent_out(i, j) = Z_indolent(i, order[j]);
    }
  }
  
  Rcpp::List out;
  out["chunk_index"] = chunk_index;
  out["n_cases"]     = n_cases;
  out["n_rows"]      = chunk_filled;
  out["Z_tau"]       = Z_tau_out;
  out["tau_accept"]  = tau_accept_out;

  if (save_Z) out["Z_indolent"] = Z_indolent_out;  
  
  std::string path = store_dir + "/" + file_prefix +
    "-latent-chunk-" + std::to_string(chunk_index) + ".rds";
  Rcout << path << std::endl;
  
  // Call back into R to write RDS
  Function saveRDS = Environment::base_env()["saveRDS"];
  saveRDS(out, Named("file") = path, Named("compress") = "xz");
  
  chunk_filled = 0;
  chunk_index += 1;
}

void MonitorLatent::finalize() {
  flush_chunk_();
}

void MonitorLatent::update(const std::vector<EPGroup_State>& groupstates,
                           const List& accept_tau,
                           const int ikept, const bool keep) {
  
  if (!keep) return;
  
  int start = 0;
  for (int k = 0; k < (int)groupstates.size(); ++k) {
    const NumericVector tmp_tau = groupstates[k].tau;
    const int n_k = tmp_tau.size();
    
    for (int j = 0; j < n_k; ++j) tau_buf[start + j] = tmp_tau[j];
    
    if (save_Z) {
      const IntegerVector tmp_Z = groupstates[k].Z;
      for (int j = 0; j < n_k; ++j) indolent_buf[start + j] = tmp_Z[j];
    }
    
    const LogicalVector tmp_acc = accept_tau[k];
    for (int j = 0; j < n_k; ++j) accept_buf[start + j] = tmp_acc[j];
    
    start += n_k;
  }
  
  for (int i = 0; i < n_cases; ++i) {
    int a = accept_buf[i];
    if (a != NA_LOGICAL) {
      n_prop[i] += 1;
      if (a) n_acc[i] += 1;
      accept_prob[i] = double(n_acc[i]) / double(n_prop[i]);
    }
  }
  
  if (save_latent == 1) {
    if (chunk_filled >= dump_every) {
      flush_chunk_();
    }
    
    Z_tau(chunk_filled, _) = tau_buf;
    tau_accept(chunk_filled, _) = accept_buf;
    if (save_Z) Z_indolent(chunk_filled, _) = indolent_buf;
    
    chunk_filled += 1;
    
    if (chunk_filled == dump_every) {
      flush_chunk_();
    }
    
  } else {
    chunk_filled = 1;
    Z_tau(0, _) = tau_buf;
    tau_accept(0, _) = accept_buf;
    if (save_Z) Z_indolent(0, _) = indolent_buf;
  }
  latest_tau = tau_buf;
  latest_ind = indolent_buf;
}
