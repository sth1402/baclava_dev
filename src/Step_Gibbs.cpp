#include <Rcpp.h>
#include "StateParameter.h"
#include "Step_Gibbs.h"
#include "Verbosity.h"

using namespace Rcpp;

// retrieve current values of parameters being updated
NumericVector Step_Gibbs::getCurrentValues() const {
  
  NumericVector value(this->np, 0.0);
  CharacterVector names(this->np);
  
  for (int i = 0; i < this->np; ++i) {
    value[i] = this->parameters[i]->value;
    names[i] = this->parameters[i]->tag;
  }
  value.names() = names;
  return value;
}

