#ifndef COMPARTMENT_H
#define COMPARTMENT_H

#include <Rcpp.h>
#include <functional>
#include <R_ext/Arith.h>
#include <variant>
#include <memory>
#include "ParameterModel.h"
#include "Distribution.h"

/* 
 * Compartment Distribution holds the full definition of a single compartment
 * level distribution, mean the Distribution used (weibull, gamma, ...) the
 * shape of that distribution, and the model assumed for the mean/rate/probability
 */

struct EPGroup_State;

class Compartment {
  
public:
  // model for the mean, rate, or probability parameter
  std::shared_ptr<ParameterModel> parameterModel;
  
  // distribution of the compartment
  std::shared_ptr<Distribution> distribution;
  
  Compartment() {};
  
  Compartment(std::shared_ptr<Distribution> distribution_,
                   std::shared_ptr<ParameterModel> parameterModel_)
    : parameterModel(std::move(parameterModel_)), distribution(std::move(distribution_)) {}
  
  Compartment(const Compartment&) = delete;
  Compartment& operator=(const Compartment&) = delete;
  
  Compartment(Compartment&&) noexcept = default;
  Compartment& operator=(Compartment&&) noexcept = default;
  
  ~Compartment() = default;
  
  // Retrieves the current values of all parameters, shape is stored last
  Rcpp::NumericVector getParameters() const;
  
  // Resets the values of all parameters, assumes shape is stored last
  void setParameters(Rcpp::NumericVector newParameterValues);
  
  // probability density function
  Rcpp::NumericVector dfunc(const Rcpp::NumericVector& x, 
                            const Rcpp::NumericVector& tau, 
                            bool log = false) const {
    const ParameterModel* pm = parameterModel.get();
    const Distribution* dist = distribution.get();
    return dist->fast_dfunc_vector(dist, x, tau, log, pm);
  }
  
  double dfunc(const double x, const double tau, bool log = false) const {
    const ParameterModel* pm = parameterModel.get();
    const Distribution* dist = distribution.get();
    return dist->fast_dfunc_scalar(dist, x, tau, log, pm);
  }
  
  // distribution function
  Rcpp::NumericVector pfunc(const Rcpp::NumericVector& x, 
                            const Rcpp::NumericVector& tau, 
                            bool lower = true, bool log = false) const {
    const ParameterModel* pm = parameterModel.get();
    const Distribution* dist = distribution.get();
    return dist->fast_pfunc_vector(dist, x, tau, lower, log, pm);
  }
  
  double pfunc(const double x, const double tau, 
               bool lower = true, bool log = false) const {
    const ParameterModel* pm = parameterModel.get();
    const Distribution* dist = distribution.get();
    return dist->fast_pfunc_scalar(dist, x, tau, lower, log, pm);
  }
  
  // quantile function
  Rcpp::NumericVector qfunc(Rcpp::NumericVector x, Rcpp::NumericVector tau, 
                      bool lower = true, bool log = false) const {
    const ParameterModel* pm = parameterModel.get();
    const Distribution* dist = distribution.get();
    return dist->qfuncf(x, tau, lower, log, pm); 
  }
  
  // backwards compatibility -- optimal scale parameter when no screening
  // data is available
  Rcpp::NumericVector compute_optimal_scale(Rcpp::NumericVector tau) const {
    const ParameterModel* pm = parameterModel.get();
    const Distribution* dist = distribution.get();
    return dist->compute_optimal_scale(tau, pm); 
  }
  
  // probability mass between a and b: P(a < X <= b)
  Rcpp::NumericVector pfunc_ab(const Rcpp::NumericVector& a, 
                               const Rcpp::NumericVector& b, 
                               const Rcpp::NumericVector& tau, 
                               bool return_log) const;
  
  // log-density adjusted for truncation to (a, b]; i.e., log of the truncated density.
  Rcpp::NumericVector dfunc_trunc(Rcpp::NumericVector x, Rcpp::NumericVector a, 
                                  Rcpp::NumericVector b, 
                                  Rcpp::NumericVector tau, bool uselog) const {
    
    const ParameterModel* pm = parameterModel.get();
    const Distribution* dist = distribution.get();
    Rcpp::NumericVector logf = dist->dfuncf(x, tau, /*log=*/true, pm);
    Rcpp::NumericVector logZ = pfunc_ab(a, b, tau, /*return_log=*/true);
    
    Rcpp::NumericVector logTrunc = logf - logZ;
    
    if (uselog) {
      return logTrunc;
    } else {
      return Rcpp::exp(logTrunc);
    }
    
  };
  
  // draw random samples from a truncated distribution
  Rcpp::NumericVector rfunc_trunc(Rcpp::NumericVector a, Rcpp::NumericVector b, 
                                  Rcpp::NumericVector tau) const {
    
    int n = a.size();
    if (b.size() != n || tau.size() != n) 
      Rcpp::stop("length mismatch: a, b, tau");
    double eps = 1e-12;
    
    const ParameterModel* pm = parameterModel.get();
    const Distribution* dist = distribution.get();
    
    Rcpp::NumericVector logFA = dist->pfuncf(a, tau, 
                                             /*lower=*/true, 
                                             /*log.p=*/true, 
                                             pm);
    Rcpp::NumericVector logFB = dist->pfuncf(b, tau, 
                                             /*lower=*/true, 
                                             /*log.p=*/true, 
                                             pm);
    
    Rcpp::NumericVector lo(n);
    Rcpp::NumericVector hi(n);
    
    for (int i = 0; i < n; ++i) {
      
      double lo_i = R_finite(logFA[i]) ? std::exp(logFA[i]) : 0.0;
      double hi_i = R_finite(logFB[i]) ? std::exp(logFB[i]) : 1.0;
      
      if (lo_i < 0.0) lo_i = 0.0;
      if (lo_i > 1.0) lo_i = 1.0;
      
      if (hi_i < lo_i) hi_i = lo_i;
      if (hi_i > 1.0) hi_i = 1.0;
      
      lo[i] = lo_i;
      hi[i] = hi_i;
    }
    
    Rcpp::LogicalVector zero = (hi <= lo + 1e-15);
    
    Rcpp::NumericVector r = Rcpp::runif(n);
    Rcpp::NumericVector u = lo + (hi - lo) * r;
    
    if (is_true(any(hi >= 1.0 | lo <= 0.0))) {
      for (int i = 0; i < lo.size(); ++i) {
        if (hi[i] >= 1.0) u[i] = std::min(u[i], 1.0 - eps);
        if (lo[i] <= 0.0) u[i] = std::max(u[i], eps);
      }
    }

    Rcpp::NumericVector q = dist->qfuncf(u, tau,
                                         /*lower=*/true,
                                         /*log.p=*/false,
                                         pm);
    
    Rcpp::LogicalVector badQ = !is_finite(q);

    if (Rcpp::is_true(any(zero))) {
      Rcpp::NumericVector zfb = Rcpp::ifelse(Rcpp::is_finite(a), a, b);
      zfb[!Rcpp::is_finite(zfb)] = NA_REAL;
      q[zero] = zfb[zero];
    }
    
    if (Rcpp::is_true(any(badQ & !zero))) {
      q[badQ & !zero] = tau[badQ & !zero];
    }
    
    return q;
  };
  
  Rcpp::NumericVector pfunc_ab_optimized(
      const EPGroup_State& gs, const Rcpp::NumericVector& tau, const double t0, 
      const Rcpp::NumericVector& wgt) const;
  
};

#endif
