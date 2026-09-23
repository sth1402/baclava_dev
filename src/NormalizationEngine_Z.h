#ifndef NORMALIZATIONENGINEZ_H
#define NORMALIZATIONENGINEZ_H

#include <Rcpp.h>
#include <RcppNumerical.h>
#include <unordered_map>
#include "Compartment.h"
#include "gauss_legendre_32.h"
#include "gauss_legendre_64.h"
#include "NegativeRateError.h"
#include "Pkg_types.h"

namespace normalization_Z {
using gl64::gauss_legendre_64;
using gl32::gauss_legendre_32;

inline void keys_by_U(const Rcpp::NumericVector& U,
                      std::vector<int>& keys, 
                      std::vector<int>& lookup) {
  std::unordered_map<int,int> seen;
  keys.reserve(U.size());
  lookup.resize(U.size());
  for (int i=0; i<U.size(); ++i) {
    int k = floor(std::round(U[i]*10.0 + 1e-12));
    auto it = seen.find(k);
    if (it == seen.end()) { 
      int idx = (int)keys.size(); 
      keys.push_back(k); 
      seen.emplace(k,idx); lookup[i]=idx;
    }
    else lookup[i] = it->second;
  }
}

struct compart {
  
  Pkg::Model model_id;
  
  bool weibull;
  bool gamma;
  bool bernoulli;
  
  // string indicating mu, rate, scale, or median
  double shape;
  Pkg::Param second_param;
  
  Rcpp::NumericVector params;
  
  double (*scale_fn)(const double, const double);
  
  compart(const Compartment& dist) {
    this->model_id = dist.parameterModel->typeid_;

    const Pkg::Dist dt = dist.distribution->typeid_;
    this->weibull = (dt == Pkg::Dist::WEIBULL);
    this->gamma = (dt == Pkg::Dist::GAMMA);
    this->bernoulli = (dt == Pkg::Dist::BERNOULLI);
    
    this->params = dist.getParameters();
    
    if (this->weibull || this->gamma) {
      this->shape = std::max(1e-6, static_cast<double>(this->params[this->params.size()-1]));
    } else {
      this->shape = NA_REAL;
    }
    
    this->second_param = dist.distribution->second_param;
    this->scale_fn = dist.distribution->scale_fn;
  }
  
  inline double scale_from(double par) const {
    return this->scale_fn(par, this->shape);
  }
};

// compute distribution mu/rate/scale value from provided parameters
inline double comp_param_fast(const compart& C, double cov) {
  switch (C.model_id) {
  case Pkg::Model::CONSTANT: return ConstantParameterModel::computeParameter(C.params, cov);
  case Pkg::Model::SIGMOID: return SigmoidParameterModel::computeParameter(C.params, cov);
  case Pkg::Model::LINEAR: return LinearParameterModel::computeParameter(C.params, cov);
  case Pkg::Model::EXP: return ExpParameterModel::computeParameter(C.params, cov);
  case Pkg::Model::STEP: return StepParameterModel::computeParameter(C.params, cov);
  case Pkg::Model::UNK: Rcpp::stop("Parameter model not recognized");
  }
  return 0.0; // unreachable
}

struct Transform {
  double U;
  double inv_beta;
  double jac_coeff;
  bool use_power;
  
  inline double map(double y) const {
    return use_power ? U * std::pow(y, inv_beta) : U * y;
  }
  
  inline double jac(double y) const {
    return use_power ? jac_coeff * std::pow(y, inv_beta - 1.0) : jac_coeff;
  }
};

struct CpIntegrand : public Numer::Func {
  
  double U;
  double t0;

  const compart* H;
  const compart* P;

  Transform tr;
  
  CpIntegrand(double U_, double t0_, 
              const compart* H_,
              const compart* P_,
              const Transform& tr_)
    :  U(U_), t0(t0_), H(H_), P(P_), tr(tr_) {}
  
  inline double dH(double x) const noexcept {
    const double par = comp_param_fast(*H, 0.0);
    const double sc = H->scale_from(par);
    
    const double shape = H->shape;
    return H->weibull ? R::dweibull(x, shape, sc, false)
                      : R::dgamma(x, shape, sc, false);
  }
  
  inline double SP(double t, double tau) const noexcept {
    const double par = comp_param_fast(*P, tau);
    const double sc = P->scale_from(par);
    
    const double shape = P->shape;
    return P->weibull ? R::pweibull(t, shape, sc, false, false)
                      : R::pgamma(t, shape, sc, false, false);
  }
  
  inline double base(double x) const noexcept {

    const double fH = dH(x);
    const double t = U - x;
    if (t < 0.0) return 0.0;
    
    const double tau = t0 + x;
    double S = SP(t, tau);
    
    return fH * S;
  }
  
  double operator()(const double& y) const override {
    if (U <= 0.0) return 0.0;
    
    const double eps = 1e-15;
    double ysafe = (y <= eps ? eps : (y >= 1.0 ? 1.0 - eps : y));
    
    const double x = tr.map(ysafe);
    const double jac = tr.jac(ysafe);
    
    return jac * base(x);
  }
};

inline Rcpp::NumericVector compute_cp_log_generic(
    Rcpp::NumericVector& AFS,
    Rcpp::IntegerVector& ind,
    double t0,
    const Compartment& Hdist,
    const Compartment& Pdist) {
  
  Rcpp::NumericVector U = AFS - t0;
  Rcpp::NumericVector zeros(U.size(), 0.0);

  Rcpp::NumericVector SH = Hdist.pfunc(U, zeros, false, false);
  if (is_true(any(!Rcpp::is_finite(SH))) || is_true(any(SH < 0.0)) ||
      is_true(any(SH > 1.0)))
    throw NegativeRateError("H pfunc returns invalid probability");

  // construct local copies of current model
  compart Pc(Pdist);
  compart Hc(Hdist);

  const double shapeH = Hc.shape;
  const bool use_power = (shapeH < 1.0 - 1e-8);

  // group by U to minimize number of integral calls
  std::vector<int> keys, lookup;
  keys_by_U(U, keys, lookup);

  Rcpp::NumericVector integ(keys.size());

  for (std::size_t i = 0; i < keys.size(); ++i) {
    double Uval = keys[i] / 10.0;

    if (!(Uval > 0.0)) {
      integ[i] = 0.0;
      continue;
    }

    Transform tr;
    tr.U = Uval;
    tr.use_power = use_power;

    if (use_power) {
      tr.inv_beta  = 1.0 / shapeH;
      tr.jac_coeff = Uval / shapeH;
    } else {
      tr.inv_beta  = 1.0;
      tr.jac_coeff = Uval;
    }

    CpIntegrand f(Uval, t0, &Hc, &Pc, tr);

    double val = use_power
      ? gauss_legendre_64([&](double y){ return f(y); })
      : gauss_legendre_32([&](double y){ return f(y); });
    if (!R_finite(val)) throw NegativeRateError("integral is not finite");

    integ[i] = val;
  }

  Rcpp::NumericVector cp(U.size());
  for (size_t i = 0; i < U.size(); ++i) cp[i] = SH[i] + integ[lookup[i]];

  for (size_t i = 0; i < cp.size(); ++i) {
    if (!(cp[i] > 0.0)) cp[i] = std::numeric_limits<double>::min();
    if (ind[i] == 1) cp[i] = 1.0;
  }
  
  return log(cp);
}
}
#endif
