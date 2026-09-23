#ifndef VERBOSITY_H
#define VERBOSITY_H

#include <Rcpp.h>
#include <iostream>
#include <string>

#pragma once
#include <iostream>
#include <string>

enum VerboseLevel {
  NONE   = 0,
  ERROR  = 1 << 0,
  INFO   = 1 << 1,
  TRACE  = 1 << 3,
  ALL    = ERROR | INFO | TRACE
};

class Verbosity {
  int mask; // which levels are active
  
public:
  Verbosity(int mask_ = NONE) : mask(mask_) {}
  
  void setMask(int m) { mask = m; }
  
  void enable(VerboseLevel lvl)   { mask |= lvl; }
  void disable(VerboseLevel lvl)  { mask &= ~lvl; }
  void clear()                    { mask = NONE; }
  bool isEnabled(VerboseLevel lvl) const { return (mask & lvl) != 0; }
  
  void log(VerboseLevel lvl, const std::string& msg) const {
    if (isEnabled(lvl)) Rcpp::Rcout << msg << std::endl;
  }
  
  void log(VerboseLevel lvl, const std::string& label, const Rcpp::NumericVector& v, int digits = 8) const {
    if (isEnabled(lvl)) {
      Rcpp::Rcout << label << to_string(v, digits) << std::endl;
    }
  }

  void log(VerboseLevel lvl, const std::string& label, const int i) const {
    if (isEnabled(lvl)) {
      Rcpp::Rcout << label << to_string_num(i) << std::endl;
    }
  }
  
  void log(VerboseLevel lvl, const std::string& label, const double d, int digits = 8) const {
    if (isEnabled(lvl)) {
      Rcpp::Rcout << label << to_string_num(d, digits) << std::endl;
    }
  }
  
  
  static std::string to_string(const Rcpp::NumericVector& v, int digits = 8) {
    std::ostringstream oss;
    oss.precision(digits);
    oss << std::fixed << "[";
    for (int i = 0; i < v.size(); ++i) {
      oss << v[i];
      if (i < v.size() - 1) oss << ", ";
    }
    oss << "]";
    return oss.str();
  }
  
  template <typename T>
  static std::string to_string_num(T x, int digits = 8) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(digits) << x;
    return oss.str();
  }
  
  
  static Verbosity& instance() {
    static Verbosity v;
    return v;
  }
  
};


#endif
