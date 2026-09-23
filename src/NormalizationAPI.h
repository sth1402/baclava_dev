#ifndef NORMALIZATIONAPI_H
#define NORMALIZATIONAPI_H

#include "NormalizationEngine.h"
#include "NormalizationEngine_Z.h"
#include "NormalizationEngine_Ztau.h"

class Compartment;

using normalization_Ztau::compute_cp_log_generic;
using normalization_Z::compute_cp_log_generic;
using normalization::compute_cp_log_generic;

// Single-component
Rcpp::NumericVector Normalization_compute_cp_log(
    Rcpp::NumericVector AFS, 
    Rcpp::IntegerVector Z,
    Rcpp::NumericVector tau,
    const Compartment& P) {
  return compute_cp_log_generic(
    AFS, Z, tau, P);
}

Rcpp::NumericVector Normalization_compute_cp_log(
    Rcpp::NumericVector AFS, 
    Rcpp::IntegerVector Z,
    double t0,
    const Compartment& H, 
    const Compartment& P) {
  return compute_cp_log_generic(
    AFS, Z, t0, H, P);
}

Rcpp::NumericVector Normalization_compute_cp_log(
    Rcpp::NumericVector AFS, 
    double t0,
    const Compartment& H, 
    const Compartment& P,
    const Compartment& psi) {
  return compute_cp_log_generic(
    AFS, t0, H, P, psi);
}


#endif
