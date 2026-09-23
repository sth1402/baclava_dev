#ifndef PRIORFUNCTIONS_H
#define PRIORFUNCTIONS_H

#include <Rcpp.h>
#include <Rmath.h>
#include <functional>
#include <map>
#include <string>
#include <vector>

class PriorFunctions {
public:
  // scalar density: f(x, params, give_log) -> double
  using PDensScalar = std::function<double(double, const std::vector<double>&, bool)>;
  
  PriorFunctions(); // default: N(0,1)
  explicit PriorFunctions(Rcpp::List R_prior); // from R list
  PriorFunctions(PDensScalar f, std::vector<double> par, std::string t);
  
  PriorFunctions(const PriorFunctions&) = default;
  PriorFunctions(PriorFunctions&&) noexcept = default;
  PriorFunctions& operator=(const PriorFunctions&) = default;
  PriorFunctions& operator=(PriorFunctions&&) noexcept = default;
  ~PriorFunctions() = default;
  
  // evaluators
  Rcpp::NumericVector log_evaluate(Rcpp::NumericVector x) const;
  double log_evaluate(double x) const;
  
  // Utilities
  inline bool is_bound() const noexcept { return static_cast<bool>(dscalar_); }
  inline void debug_dump() const;
  
  // type and parameter settings
  const std::string& type()   const noexcept { return type_; }
  const std::vector<double>& params() const noexcept { return params_; }
  
  inline std::string getType() const {return type_;}
  
private:
  PDensScalar dscalar_;
  std::vector<double> params_;
  std::string type_;
  
  // Helpers
  static PDensScalar bind_scalar_for_type(const std::string& dist);
  static int required_param_min(const std::string& dist);
  static void validate_params_or_throw(const std::string& dist,
                                       const std::vector<double>& p);
};

#endif // PRIORFUNCTIONS_H
