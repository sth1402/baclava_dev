#ifndef COVADAPTER_H
#define COVADAPTER_H

#include <Rcpp.h>
#include <vector>
#include <memory>
#include <limits>
#include "Verbosity.h"

/*
 * CovAdapter
 */

struct CovAdapter {
  
  CovAdapter() {};
  
  explicit CovAdapter(const Rcpp::List& adapt_list)  {
    
    Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),
                  "[covadapt] initializing from list");
    
    if (!adapt_list.containsElementNamed("meta_for"))
      Rcpp::stop("CovAdapter Constructor: provided Adapt list invalid");
    
    if (Rcpp::as<std::string>(adapt_list["meta_for"]) != "Adapt")
      Rcpp::stop("CovAdapter Constructor: provided list is not of type Adapt");
    
    if (!adapt_list.containsElementNamed("ridge") ||
        !adapt_list.containsElementNamed("learn"))
        Rcpp::stop("adapt_list is not properly formatted");
    
    try{
      this->ridge = Rcpp::as<double>(adapt_list["ridge"]);

      std::string learning = Rcpp::as<std::string>(adapt_list["learn"]);
      this->learn_L = learning == "sDL";
      this->learn_D = learning == "sD" || learning == "sDL";
    } catch (const std::exception& e) {
      Rcpp::stop("CovAdapter constructor: missing or malformed parameter list: %s", e.what());
    }
    if (!std::isfinite(ridge) || ridge <= 0.0) Rcpp::stop("ridge must be > 0");
    if (!std::isfinite(ridge_max) || ridge_max < ridge) Rcpp::stop("ridge_max must be >= ridge");
    
    if (this->learn_L && !this->learn_D) Rcpp::stop("must learn D or D and L");

    Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO),
                  "[covadapt] ridge=" + std::to_string(this->ridge)
                  + " ridge_max=" + std::to_string(this->ridge_max)
                  + " learn L=" + std::to_string(this->learn_L)
                  + " learn D=" + std::to_string(this->learn_D));
  }
  
  CovAdapter(const CovAdapter&) = default;
  CovAdapter& operator=(const CovAdapter&) = default;
  CovAdapter(CovAdapter&&) = default;
  CovAdapter& operator=(CovAdapter&&) = default;
  
  void update(const Rcpp::NumericVector& x, 
              Rcpp::NumericMatrix& L, 
              Rcpp::NumericVector& D, 
              const int m, const bool run_covD, const bool run_covL) {
    
    if (this->frozen) return;

    Verbosity::instance().log(TRACE,
                  "[covadapt] " + std::to_string(m));
    
    welford_update(x, L, D, run_covD, run_covL);
    
  }
  
  bool isFrozen() { return this->frozen; }
  void freeze() { this->frozen = true; }
  void request_freeze() { this->frozen = true; }
  
private:
  
  // at what iteration to end updates
  // is covariance matrix frozen
  bool frozen = false;
  bool learn_L = true;
  bool learn_D = true;
  
  // running stats on WORKING scale
  Rcpp::NumericMatrix M2; // d x d sum of outer products for covariance
  Rcpp::NumericVector mean; // size d
  
  int n_seen = 0;
  int chol_start = 200;
  int chol_count = 50;

  // adaptation bookkeeping
  double ridge = 1e-8; // penalty of Welford update
  double ridge_max = 1e-2;
  
  // ensure internal arrays are sized to d and zero-initialized where needed
  void ensure_dims(int d) {
    if (this->M2.nrow() != d || this->M2.ncol() != d) this->M2 = Rcpp::NumericMatrix(d, d);
    if (this->mean.size() != d) this->mean = Rcpp::NumericVector(d);
  }
  
  // compute the lower-triangular Cholesky factor L of a symmetric
  // positive-definite matrix A such that A = L * L^T.
  //
  // On success, returns true and writes the factor to Lout.
  // On failure, returns false and leaves Lout in its current state.
  bool chol_lower(const Rcpp::NumericMatrix& A, Rcpp::NumericMatrix& Lout) {
    const int d = A.nrow();
    if (A.ncol() != d) return false;
    
    Lout = Rcpp::NumericMatrix(d, d);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j <= i; ++j) {
        double s = A(i, j);
        for (int k = 0; k < j; ++k) s -= Lout(i, k) * Lout(j, k);
        if (i == j) {
          if (!R_finite(s) || s <= 0.0) return false;
          Lout(i, j) = std::sqrt(s);
        } else {
          if (Lout(j, j) == 0.0) return false;
          Lout(i, j) = s / Lout(j, j);
        }
      }
    }
    return true;
  }
  
  Rcpp::NumericMatrix cov() const {
    const int d = this->M2.nrow();
    Rcpp::NumericMatrix C(d, d);
    
    if (n_seen < 2) {
      for (int i = 0; i < d; ++i) C(i,i) = 1e-6;
      return C;
    }
    
    const double denom = static_cast<double>(n_seen - 1);
    for (int i = 0; i < d; ++i)
      for (int j = 0; j < d; ++j)
        C(i,j) = this->M2(i,j) / denom;
    
    for (int i = 0; i < d; ++i)
      for (int j = i + 1; j < d; ++j) {
        const double a = 0.5 * (C(i,j) + C(j,i));
        C(i,j) = C(j,i) = a;
      }
      
    return C;
  }
 
 /*
  * Attempt to generates a lower-triangle Cholesky factor (L) from covariance matrix.
  * Returns true of cholesky decomposition successful; false otherwise
  */ 
  void refresh_chol(Rcpp::NumericMatrix& L, Rcpp::NumericVector& D, const bool run_covD, const bool run_covL) {
    
    Rcpp::NumericMatrix Lnew;
    
    Rcpp::NumericMatrix C = cov();
    const int d = C.nrow();
    
    // clamp diag to avoid zeros
    if (this->learn_D && run_covD) {
      double mean_diag = 0.0;
      for (int i = 0; i < d; ++i) mean_diag += C(i,i);
      mean_diag /= d;
      mean_diag = std::max(mean_diag, 0.0);
      const double min_var = std::max(1e-12, 1e-6 * mean_diag); // allow eps to be input later
      for (int i = 0; i < d; ++i) {
        double v = C(i, i);
        if (!R_finite(v) || v < min_var) v = min_var;
        D[i] = std::sqrt(v);
      }
    }
    
    if (!this->learn_L) return;
    
    // Build correlation R
    Rcpp::NumericMatrix R(d, d);
    for (int i = 0; i < d; ++i) {
      R(i,i) = 1.0;
      for (int j = i+1; j < d; ++j) {
        double denom = D[i] * D[j];
        double rij = (denom > 0.0) ? C(i,j) / denom : 0.0;
        if (!R_finite(rij)) rij = 0.0;
        rij = std::max(-0.999, std::min(0.999, rij));
        R(i,j) = R(j,i) = rij;
      }
    }
    
    if (chol_lower(R, Lnew)) {
      if (run_covL) L = Lnew;
      std::string msg = "[corrBlock] chol no ridge=";
      for (int u = 0; u < d; ++u) msg = msg + " " + std::to_string(L(u,u));
      Verbosity::instance().log(TRACE, msg);
      return;
    }
    
    
    double jitter = this->ridge;
    while (jitter <= this->ridge_max) {
      Rcpp::NumericMatrix A = Rcpp::clone(R);
      for (int i = 0; i < d; ++i) A(i,i) += jitter;
      
      if (chol_lower(A, Lnew)) {
        if (run_covL) L = Lnew;
        std::string msg = "[corrBlock] chol ridge=";
          for (int u = 0; u < d; ++u) msg = msg + " " + std::to_string(L(u,u));
          Verbosity::instance().log(TRACE, msg);
          return;
      }
      jitter *= 10.0;
    }
    
    for (int i = 0; i < d; ++i)
      for (int j = 0; j < d; ++j)
        L(i,j) = (i == j) ? 1.0 : 0.0;
    
    for (int i = 0; i < d; ++i) 
      if (!R_finite(D[i])) { 
        this->frozen = true; 
        std::string msg = "[corrBlock] chol frozen due to pathology";
        Rcpp::warning(msg);
        Verbosity::instance().log(TRACE, msg);
        return; 
      }
  }
  
  void welford_update(const Rcpp::NumericVector& x,
                      Rcpp::NumericMatrix& L,
                      Rcpp::NumericVector& D,
                      const bool run_covD, const bool run_covL) {
    
    const int d = x.size();
    ensure_dims(d);

    if (n_seen == 0) {
      mean = clone(x);
      for (int i = 0; i < d; ++i) 
        for (int j = 0; j < d; ++j) {
          this->M2(i, j) = 0.0;
        }
      n_seen = 1;
      return;
    }

    n_seen++;
    
    // running mean
    Rcpp::NumericVector delta = x - mean;
    this->mean = mean + delta / n_seen;
    Rcpp::NumericVector delta2 = x - this->mean;
    // running Sigma
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        if (!this->learn_L && i != j) continue;
        this->M2(i, j) += delta[i] * delta2[j];
      }
    }
    
    if (n_seen >= chol_start && ((n_seen % chol_count) == 0)) {
      refresh_chol(L, D, run_covD, run_covL);
      return;
    }

  }
};

#endif // COVADAPTER_H
