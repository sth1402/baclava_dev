#ifndef UPDATES_H
#define UPDATES_H

#include <Rcpp.h>

class Model;
class Step_MH;
class Step_Gibbs;
struct EPGroup_State;

class Updates {
public:
  static void H(std::vector<EPGroup_State>& groupstates, Model* model, int m);
  static void MH(std::vector<EPGroup_State>& groupstates, Model* model, Step_MH& group, int m);
  
  static void P(std::vector<EPGroup_State>& groupstates, Model* model, int m);
  
  // joint update of psi and Z_indolent
  static void psiZ(std::vector<EPGroup_State>& groupstates, Model* model, Step_MH& group, int m);
  
  static void beta(std::vector<EPGroup_State>& groupstates,
                   Model* model,
                   Step_Gibbs& group);

  static void tau(std::vector<EPGroup_State>& groupstates, 
                  Rcpp::List& accepted, 
                  Model* model);
};
#endif
