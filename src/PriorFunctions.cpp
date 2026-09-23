#include "PriorFunctions.h"
#include "Verbosity.h"
#include <stdexcept>

// prior parameter requirements
int PriorFunctions::required_param_min(const std::string& d) {
  if (d == "no_prior")  return 0;
  if (d == "dunif")     return 2;              // min, max
  if (d == "dnorm")     return 2;              // mean, sd
  if (d == "dhalfnorm") return 1;              // sd
  if (d == "dlnorm")    return 2;              // meanlog, sdlog
  if (d == "dgamma")    return 2;              // shape, scale   (NOTE: scale, not rate)
  if (d == "dinvgamma") return 2;              // shape, scale   (NOTE: scale, not rate)
  if (d == "dbeta")     return 2;              // a, b
  if (d == "dchisq")    return 1;              // df
  if (d == "dnchisq")   return 2;              // df, ncp
  if (d == "dt")        return 1;              // df
  if (d == "dnt")       return 2;              // df, ncp
  if (d == "df")        return 2;              // df1, df2
  if (d == "dnf")       return 3;              // df1, df2, ncp
  if (d == "dcauchy")   return 2;              // loc, scale
  if (d == "dexp")      return 1;              // rate
  if (d == "dlogis")    return 2;              // loc, scale
  if (d == "dweibull")  return 2;              // shape, scale
  return -1; // unknown
}

void PriorFunctions::validate_params_or_throw(const std::string& dist,
                                              const std::vector<double>& p) {
  
  // ensure p is of appropriate length
  const int mn = required_param_min(dist);
  if (mn < 0) Rcpp::stop("unsupported prior distribution '%s'.", dist.c_str());
  
  if ((int)p.size() < mn)
    Rcpp::stop("prior '%s' expects at least %d parameter(s); got %d.",
               dist.c_str(), mn, (int)p.size());
}

// binding table
PriorFunctions::PDensScalar PriorFunctions::bind_scalar_for_type(const std::string& dist) {
  
  if (dist == "no_prior") {
    return [](double, const std::vector<double>&, bool give_log) { return give_log ? 0.0 : 1.0; };
  }
  
  if (dist == "dunif") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dunif(x, a[0], a[1], l); };
  
  if (dist == "dnorm") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dnorm4(x, a[0], a[1], l); };
  
  if (dist == "dhalfnorm") return [=](double x, const std::vector<double>& a, bool l) {
    if (x < 0.0) return l ? R_NegInf : 0.0;
    double dens = R::dnorm4(x, 0.0, a[0], l);
    return l ? dens + std::log(2.0) : 2.0 * dens;
  };  
  
  if (dist == "dlnorm") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dlnorm(x, a[0], a[1], l); };
  
  if (dist == "dgamma") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dgamma(x, a[0], a[1], l); };
  
  if (dist == "dinvgamma") return [=](double x, const std::vector<double>& a, bool l) { 
    if (x <= 0.0) return l ? R_NegInf : 0.0;
    double logdens =
      - a[0] * log(a[1])
      - R::lgammafn(a[0])
      - (a[0] + 1.0) * log(x)
      - 1.0 / (a[1] * x);
    return l ? logdens : std::exp(logdens);
  };
  if (dist == "dbeta") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dbeta(x, a[0], a[1], l); };
  
  if (dist == "dchisq") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dchisq(x, a[0], l); };
  
  if (dist == "dnchisq") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dnchisq(x, a[0], a[1], l); };
  
  if (dist == "dt") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dt (x, a[0], l); };
  
  if (dist == "dnt") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dnt(x, a[0], a[1], l); };
  
  if (dist == "df") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::df (x, a[0], a[1], l); };
  
  if (dist == "dnf") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dnf(x, a[0], a[1], a[2], l); };
  
  if (dist == "dcauchy") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dcauchy(x, a[0], a[1], l); };
  
  if (dist == "dexp") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dexp(x, a[0], l); };
  
  if (dist == "dlogis") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dlogis(x, a[0], a[1], l); };
  
  if (dist == "dweibull") return [=](double x, const std::vector<double>& a, bool l) { 
    return R::dweibull(x, a[0], a[1], l); };
  
  // unknown
  return PDensScalar{};
}

// default constructor
PriorFunctions::PriorFunctions()
  : dscalar_(bind_scalar_for_type("dnorm")),
    params_({0.0, 1.0}),
    type_("dnorm") {}

// constructor from List defined in R's Prior class
PriorFunctions::PriorFunctions(Rcpp::List R_prior) : PriorFunctions() {
  if (R_prior.size() == 0) return;
  
  if (!R_prior.containsElementNamed("meta_for") ||
      Rcpp::as<std::string>(R_prior["meta_for"]) != "Prior")
    Rcpp::stop("PriorFunctions: provided list is not of type 'Prior'.");
  
  if (!R_prior.containsElementNamed("dist"))
    Rcpp::stop("PriorFunctions: missing field 'dist'.");
  
  const std::string dist_type = Rcpp::as<std::string>(R_prior["dist"]);
  
  if (dist_type == "no_prior") {
    type_ = "no_prior";
    params_.clear();
    dscalar_ = bind_scalar_for_type(type_);
    return;
  }
  
  if (!R_prior.containsElementNamed("params"))
    Rcpp::stop("PriorFunctions: missing field 'params' for dist '%s'.", dist_type.c_str());
  
  type_ = dist_type;
  
  Rcpp::NumericVector par = R_prior["params"];
  params_.assign(par.begin(), par.end());
  
  validate_params_or_throw(type_, params_);
  
  dscalar_ = bind_scalar_for_type(type_);
  if (!dscalar_)
    Rcpp::stop("PriorFunctions: unsupported dist '%s'.", type_.c_str());
}

PriorFunctions::PriorFunctions(PDensScalar f, std::vector<double> par, std::string t)
  : dscalar_(std::move(f)), params_(std::move(par)), type_(std::move(t)) {
  validate_params_or_throw(type_, params_);
  if (!dscalar_) Rcpp::stop("PriorFunctions: empty function pointer for '%s'.", type_.c_str());
}

Rcpp::NumericVector PriorFunctions::log_evaluate(Rcpp::NumericVector x) const {
  Rcpp::NumericVector out(x.size());
  for (R_xlen_t i = 0; i < x.size(); ++i)
    out[i] = dscalar_(x[i], params_, /*give_log=*/true);
  return out;
}

double PriorFunctions::log_evaluate(double x) const {
  double out = dscalar_(x, params_, /*give_log=*/true);
  return out;
}

void PriorFunctions::debug_dump() const {
  std::string msg = "[Prior] type=" + type_ 
    + " params_len=" + std::to_string(params_.size())
    + " bound=" + std::to_string((bool)dscalar_)  + " params=";
  if (params_.empty()) msg = msg + "(none)";
  else {
    msg = msg + "{";
    for (int i=0;i<params_.size();++i) {
      if (i) msg = msg + ",";
      msg = msg + std::to_string(params_[i]);
    }
    msg = msg + "}";
  }
  Verbosity::instance().log(static_cast<VerboseLevel>(TRACE | INFO), msg);
}
