#ifndef STATEPARAMETER_H
#define STATEPARAMETER_H

#include <Rcpp.h>
#include <cmath>
#include <memory>
#include "PriorFunctions.h"
#include "ProposalScale.h"
#include "NegativeRateError.h"
#include "Pkg_types.h"

/*
 * Primary data structure of every model parameter. 
 * Includes a name, the distribution to which it belongs, its current value,
 * if it is fixed or updated, the evaluation of its log-prior, how to 
 * handle proposals (transformation/reflections), and any Jacobian contributions
 */

class StateParameter {
public:
  
  // user specified name for the parameter
  std::string tag = "unk";
  
  // user specified name for the distribution to which the parameter pertains
  std::string from_dist_parameter = "unk";
  
  Pkg::Components component = Pkg::Components::UNK;
  
  // current value of the parameter
  double value = 0.0;
  
  // if true, do not allow parameter value to change
  bool fixed = false;
  
  std::string prop_func = "rnorm";

  // update value
  void update(double new_value) {
    if (fixed) Rcpp::stop("cannot update a fixed parameter");
    
    if ((ps->tag != "id") && (new_value <= minv || new_value >= maxv))
      throw NegativeRateError("value at or beyond boundary of support");
    if ((ps->tag == "id") && (new_value < minv || new_value > maxv)) 
      throw NegativeRateError("value out of bounds");
    
    value = new_value;
    for (auto it = mirrors.begin(); it != mirrors.end(); ) {
      if (auto m = it->lock()) { m->value = new_value; ++it; }
      else { it = mirrors.erase(it); }
    }
  }
  
  // evaluate log-prior at x
  double log_prior(double x) const { 
    return prior.log_evaluate(x); 
  }
  
  Rcpp::NumericVector log_prior(Rcpp::NumericVector x) const { 
    return prior.log_evaluate(x); 
  }
  
  // retrieve current value
  double get_value() const { return value; }
  
  // retrieve current value with its tag
  Rcpp::List get_tag_value() const { 
    return Rcpp::List::create(Rcpp::Named("tag") = tag, Rcpp::Named("value") = value);
  }
  
  /*
   * when one subset allows for fitting the shape, but the other subset shape
   * is fixed, we mirror the shape in both subsets.
   * This is for Y/O exploration and probably should not be accessible in
   * the package
   */
  void add_mirror(const std::shared_ptr<StateParameter>& m) {
    mirrors.emplace_back(m);
    m->fixed = true;
  }
  
  // default constructor
  StateParameter() : 
    tag("unk"), from_dist_parameter("unk"), component(Pkg::Components::UNK),
    value(0.0), fixed(true), prop_func("rnorm"), 
    minv(-INFINITY), maxv(INFINITY), prior() {
    make_scale("id");
  }
  
  std::shared_ptr<StateParameter> clone() const {
    return std::make_shared<StateParameter>(*this);
  }
  
  /*
   * Param$toList = function() {
   *   list("meta_for" = "Param",
   *        "tag" = self$tag, 
   *.       "init" = self$init,
   *        "minv" = self$minv,
   *        "maxv" = self$maxv,
   *        "prop_scale" = self$prop_scale, 
   *        "prop_func" = self$prop_func,
   *        "prior" = self$prior$toList(),
   *        "fixed" = self$fixed)
   },
   */
  StateParameter(Rcpp::List parameter, Pkg::Components component_) {
    
    if (!parameter.containsElementNamed("meta_for"))
      Rcpp::stop("StateParameter constructor: provided parameter list invalid");
    
    if (Rcpp::as<std::string>(parameter["meta_for"]) != "Param")
      Rcpp::stop("StateParameter constructor: provided list is not of type Param");

    if (!parameter.containsElementNamed("init") ||
        !parameter.containsElementNamed("tag") ||
        !parameter.containsElementNamed("minv") ||
        !parameter.containsElementNamed("maxv") ||
        !parameter.containsElementNamed("prop_scale") ||
        !parameter.containsElementNamed("prop_func") ||
        !parameter.containsElementNamed("prior") ||
        !parameter.containsElementNamed("fixed"))
        Rcpp::stop("parameter is not properly formatted");
    
    try {
      
      this->tag = Rcpp::as<std::string>(parameter["tag"]);
      this->from_dist_parameter = "";
      this->value = Rcpp::as<double>(parameter["init"]);
      this->minv = Rcpp::as<double>(parameter["minv"]);
      this->maxv = Rcpp::as<double>(parameter["maxv"]);
      this->fixed = Rcpp::as<bool>(parameter["fixed"]);
      this->prior = PriorFunctions(Rcpp::as<Rcpp::List>(parameter["prior"]));
      std::string scale = Rcpp::as<std::string>(parameter["prop_scale"]);
      make_scale(scale);
      this->prop_func = Rcpp::as<std::string>(parameter["prop_func"]);
      
      this->component = component_;

    } catch (const std::exception& e) {
      Rcpp::stop("StateParameter constructor: missing or malformed parameter list: %s", e.what());
    }
  }
  
  StateParameter(const StateParameter& other)
    : tag(other.tag), from_dist_parameter(other.from_dist_parameter),
      component(other.component),
      value(other.value), prop_func(other.prop_func), 
      minv(other.minv), maxv(other.maxv), prior(other.prior) {
    ps = other.ps->clone();
  }
  
  StateParameter(StateParameter&& other) noexcept
    : tag(std::move(other.tag)),
      from_dist_parameter(std::move(other.from_dist_parameter)),
      component(std::move(other.component)),
      value(std::move(other.value)), 
      fixed(std::move(other.fixed)),
      prop_func(std::move(other.prop_func)),
      minv(std::move(other.minv)),
      maxv(std::move(other.maxv)),
      prior(std::move(other.prior)),
      ps(std::move(other.ps)) {}
  
  // Copy assignment operator
  StateParameter& operator=(const StateParameter& other) {
    if (this == &other) return *this;
    
    this->tag = other.tag;
    this->from_dist_parameter = other.from_dist_parameter;
    this->component = other.component;
    this->value = other.value;
    this->prior = other.prior;
    this->ps = other.ps->clone();
    this->fixed = other.fixed;
    this->minv = other.minv;
    this->maxv = other.maxv;
    this->prop_func = other.prop_func;

    return *this;
  }
  
  StateParameter& operator=(StateParameter&& other) noexcept {
    if (this == &other) return *this;
    
    this->tag = std::move(other.tag);
    this->from_dist_parameter = std::move(other.from_dist_parameter);
    this->component = std::move(other.component);
    this->value = std::move(other.value);
    this->prior = std::move(other.prior);
    this->ps = std::move(other.ps);
    this->fixed = std::move(other.fixed);
    this->minv = std::move(other.minv);
    this->maxv = std::move(other.maxv);
    this->prop_func = std::move(other.prop_func);

    return *this;
  }
  
  ~StateParameter() = default;
  
  double to_working(double value_nat) const { return this->ps->apply_scale(value_nat); }
  
  double from_working(double value_work) const { return this->ps->invert_scale(value_work); }

  /*
   * log of the absolute Jacobian determinant of the theta -> w transformation
   * theta is the parameter on the natural scale, w is the working (proposal)
   * scale
   */
  inline double log_abs_dthetadw(double theta) const {
    return this->ps->log_abs_dxdw(theta);
  }
  
  PriorFunctions getPrior() const { return prior; }
  
  double apply_boundaries(double& x) {
    if (ps->tag == "logit" || 
        ps->tag == "log_pos" ||
        ps->tag == "log_shifted" || 
        ps->tag == "log_flipped")
      return x;    
    
    double lo = this->minv;
    double hi = this->maxv;
    
    int cnt = 0;
    while (true) {
      if (cnt > 100) throw NegativeRateError("reflection failed");
      cnt++;
      if (std::isfinite(lo) && x < lo) { x = 2.0*lo - x; continue; }
      if (std::isfinite(hi) && x > hi) { x = 2.0*hi - x; continue; }
      break;
    }
    return x;
  }
  
private:
  
  double minv = -INFINITY;
  double maxv =  INFINITY;
  
  // prior for parameter update decisions
  PriorFunctions prior;
  
  // ProposalScale object to handle all of the w <-> theta logic
  std::unique_ptr<ProposalScale> ps;

  // allowing for subsets to share a shape parameter
  std::vector<std::weak_ptr<StateParameter>> mirrors;
  
  void make_scale(const std::string& s) {
    if (s == "natural" || s == "id") {
      this->ps = std::make_unique<NaturalScale>();
    } else if (s == "log_pos") {
      this->ps =  std::make_unique<PosLogScale>(); 
    } else if (s == "log_shifted") {
      this->ps =  std::make_unique<ShiftedLogScale>(this->minv); 
    } else if (s == "log_flipped") {
      this->ps =  std::make_unique<FlippedLogScale>(this->maxv); 
    } else if (s == "logit") {
      this->ps =  std::make_unique<LogitScale>(this->minv, this->maxv); 
    } else {
      throw std::runtime_error("unknown scale: " + s);
    }
  }
 
};

#endif
