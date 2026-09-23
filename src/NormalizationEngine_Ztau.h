#ifndef NORMALIZATIONENGINEZTAU_H
#define NORMALIZATIONENGINEZTAU_H

#include <Rcpp.h>
#include <RcppNumerical.h>
#include <unordered_map>
#include "Compartment.h"
#include "NegativeRateError.h"
#include "Pkg_types.h"

namespace normalization_Ztau {


inline Rcpp::NumericVector compute_cp_log_generic(
    Rcpp::NumericVector& AFS,
    Rcpp::IntegerVector& ind,
    Rcpp::NumericVector& tau,
    const Compartment& Pdist) {

  Rcpp::NumericVector cp(AFS.size(), 1.0);
  
  for (int i = 0; i < AFS.size(); ++i) {
    if (ind[i] == 0 && tau[i] < AFS[i]) {
      cp[i] = Pdist.pfunc(AFS[i] - tau[i], 0.0, false, false);
    }
  }

  return log(cp);
}
}
#endif
