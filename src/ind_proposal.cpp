#include <Rcpp.h>
#include <cmath>
#include "Model.h"
#include "ind_proposal.h"
#include "EPGroup_State.h"
#include "Pkg_types.h"

using namespace Rcpp;

namespace {

inline double clamp01(const double p, double eps = 1e-12) {
  double out = p;
  if (!R_finite(p)) out = 0.5;       // defensive, should rarely trigger
  out = std::max(eps, std::min(1.0 - eps, out));
  return out;
}

inline double lse2(double a, double b) {
  if (a > b) return a + std::log1p(std::exp(b - a));
  else return b + std::log1p(std::exp(a - b));
}

} //end namespace

// Compute the probability of being indolent for each individual,
// given their current tau values and the model parameters.
//
// This is done by computing the likelihood under both possible indolence states:
//    L_0 = likelihood assuming progressive tumor (ind = 0)
//    L_1 = likelihood assuming indolent tumor    (ind = 1)
// and then returning the probability:
//    P(indolent = 1 | tau, data, model) = L_1 / (L_0 + L_1)
NumericVector compute_prob_indolent(const EPGroup_State& gs, 
                                    const NumericVector& tau, 
                                    const Model* model) {
  
  uint8_t impact = 0;
  impact |= static_cast<uint8_t>(Pkg::Components::Z);
  
  NumericVector result(gs.n, 0.0);
  if (gs.endpoint_type == Pkg::EP::CLINICAL) return result;

  IntegerVector ind(gs.n, 0);
  NumericVector L_0 = gs.log_likelihood(model, ind, impact);

  ind.fill(1);
  NumericVector L_1 = gs.log_likelihood(model, ind, impact);

  for (int i = 0; i < gs.n; ++i) {
    const double denom = lse2(L_0[i], L_1[i]);
    result[i] = std::exp(L_1[i] - denom);
    result[i] = clamp01(result[i], 1e-12);
  }

  return result;
}

IntegerVector rprop_indolent(const Pkg::EP endpoint_type, 
                             const NumericVector& prob_indolent) {
  
  int n = prob_indolent.size();
  IntegerVector indolent(n, NA_INTEGER);
  
  if (endpoint_type == Pkg::EP::CLINICAL) {
    // Clinical cases
    // indolence is always false for those with a clinical diagnosis
    indolent.fill(0);
  } else {
    // Censored and preclinical cases
    // indolence is a random sample from binomial using the current
    // estimated probability of indolence
    for (int i = 0; i < n; ++i) {
      indolent[i] = R::rbinom(1, prob_indolent[i]);
    }
  }
  
  return indolent;
}

// vector of dimension equal to the number of participants
NumericVector dlog_prop_indolent(const Pkg::EP endpoint_type, 
                                 const NumericVector& prob_indolent,
                                 const IntegerVector& indolent) {
  
  int n = prob_indolent.size();
  
  NumericVector result(n, 0.0);
  if (endpoint_type == Pkg::EP::CLINICAL) {
    return result;
  } else {
    for (int i = 0; i < n; ++i) {
      result[i] = indolent[i] == 1 ? std::log(prob_indolent[i]) : std::log1p(-prob_indolent[i]);
    }
    return result;
  }
}
