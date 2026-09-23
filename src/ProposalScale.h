#ifndef PROPOSALSCALE_H
#define PROPOSALSCALE_H

#include <cmath>
#include <limits>
#include <memory>
#include "NegativeRateError.h"

class ProposalScale {
public:
  
  std::string tag;
  
  virtual ~ProposalScale() {}
  
  virtual double apply_scale(double x) const {
    return x;
  };
  
  virtual double invert_scale(double w) const {
    return w;
  }
  
  virtual double log_abs_dxdw(double x) const {
    return 0.0;
  }
  
  virtual std::unique_ptr<ProposalScale> clone() const = 0;
  
};

class NaturalScale : public ProposalScale {
public:
  
  NaturalScale() {
    this->tag = "id";
  }
  
  double apply_scale(double x) const override { return x; }
  
  double invert_scale(double w) const override { return w; }
  
  double log_abs_dxdw(double x) const override { return 0.0; }
  
  std::unique_ptr<ProposalScale> clone() const override {
    return std::make_unique<NaturalScale>(*this);
  }
};

// Support is positive axis theta \in (0, Inf)
class PosLogScale : public ProposalScale {
public:
  
  PosLogScale() { this->tag = "log_pos"; };
  
  double apply_scale(double x) const override {
    if (!std::isfinite(x)) throw NegativeRateError("infinite parameter on the natural scale");
    if (x <= 0.0) throw NegativeRateError("theta <= 0 bound in PosLogScale");
    return std::log(x);
  };
  
  double invert_scale(double w) const override {
    if (!std::isfinite(w)) 
      throw NegativeRateError("infinite parameter on the working scale");
    w = std::max(w, -700.0); 
    return std::exp(w); 
  }
  
  double log_abs_dxdw(double x) const override {
    if (x <= 0) Rcpp::stop("log-scale update requires theta > 0");
    return std::log(x);
  }
  
  std::unique_ptr<ProposalScale> clone() const override {
    return std::make_unique<PosLogScale>(*this);
  }
  
};

// Support is theta \in (a, Inf)
class ShiftedLogScale : public ProposalScale {
public:
  
  ShiftedLogScale(double a_) : a(a_) { 
    this->tag = "log_shifted";
    if (!std::isfinite(a_))
      throw std::runtime_error("ShiftedLogScale requires finite lower bound");
  };
  
  double apply_scale(double x) const override {
    if (!std::isfinite(x))
      throw NegativeRateError("non-finite theta in ShiftedLogScale");
    
    if (x <= this->a)
      throw NegativeRateError("theta <= lower bound in ShiftedLogScale");
    
    return std::log(x - this->a);
  };
  
  double invert_scale(double w) const override {
    if (!std::isfinite(w))
      throw NegativeRateError("non-finite w in ShiftedLogScale");
    w = std::max(w, -700.0);
    return this->a + std::exp(w);
  }
  
  double log_abs_dxdw(double x) const override {
    if (x - this->a <= 0) Rcpp::stop("log-scale update requires theta > 0");
    return std::log(x - this->a);
  }
  
  std::unique_ptr<ProposalScale> clone() const override {
    return std::make_unique<ShiftedLogScale>(*this);
  }
  
private:
  
  double a = 0.0;
  
};

class FlippedLogScale : public ProposalScale {
public:
  explicit FlippedLogScale(double b_) : b(b_) {
    this->tag = "log_flipped";
    if (!std::isfinite(b))
      throw std::runtime_error("FlippedLogScale requires finite upper bound");
  }
  
  double apply_scale(double x) const override {
    if (!std::isfinite(x))
      throw NegativeRateError("non-finite theta in FlippedLogScale");
    
    if (x >= b)
      throw NegativeRateError("theta >= upper bound in FlippedLogScale");
    
    return std::log(b - x);
  }
  
  double invert_scale(double w) const override {
    if (!std::isfinite(w))
      throw NegativeRateError("non-finite w in FlippedLogScale");
    w = std::max(w, -700.0);
    return b - std::exp(w);
  }
  
  double log_abs_dxdw(double x) const override {
    if (x >= b)
      throw NegativeRateError("theta >= upper bound in FlippedLogScale");
    return std::log(b - x);
  }
  
  std::unique_ptr<ProposalScale> clone() const override {
    return std::make_unique<FlippedLogScale>(*this);
  }
  
private:
  double b;
};


class LogitScale : public ProposalScale {
public:
  
  LogitScale() : LogitScale(0.0, 1.0) {}
  
  LogitScale(double a_, double b_) : a(a_), b(b_) {
    this->tag = "logit";
    
    if (!std::isfinite(a_) || !std::isfinite(b_))
      throw std::runtime_error("logit scale requires finite [a,b] bounds");
    
    if (b_ <= a_)
      throw std::runtime_error("logit scale requires b > a");
  };
  
  double apply_scale(double x) const override {

    if (!std::isfinite(x)) throw NegativeRateError("infinite or NaN parameter on the logit scale");
    if (x <= this->a || x >= this->b) throw NegativeRateError("infinite or NaN parameter on the logit scale");
    
    double t = (x - this->a) / (this->b - this->a);
    return std::log(t) - std::log1p(-t);
  };
  
  double invert_scale(double w) const override {
    if (!std::isfinite(w)) throw NegativeRateError("infinite or NaN parameter on the working scale");
    w = std::max(w, -700.0);
    w = std::min(w,  700.0);
    return this->a + (this->b - this->a) / (1.0 + std::exp(-w));
  }
  
  double log_abs_dxdw(double x) const override {
    
    if (x <= this->a || x >= this->b)
      return -std::numeric_limits<double>::infinity();
    
    const double eps = 1e-12;
    double z = (x - this->a) / (b - this->a);
    z = std::min(std::max(z, eps), 1.0 - eps);
    
    return std::log(this->b - this->a) + std::log(z) + std::log(1.0 - z);
  }
  
  std::unique_ptr<ProposalScale> clone() const override {
    return std::make_unique<LogitScale>(*this);
  }
  
private:
  double a = 0.0;
  double b = 1.0;
};

#endif
