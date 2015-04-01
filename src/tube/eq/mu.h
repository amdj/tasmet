// mu.h
//
// Author: J.A. de Jong 
//
// Description:
// Equation describing the relation between momentum flow at the
// vertex and the mass flow at the cell walls and density at the vertex
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef MU_H
#define MU_H
#include "tubeequation.h"

namespace tasystem {
  class JacRow;
} // namespace tasystem

namespace tube {
  SPOILNAMESPACE  
  class Cell;
  class Tube;
  class WeightFactors;

  class Mu: public Equation{
  public:
    Mu(const Cell& v):Equation(v){}
    vd error() const;

    tasystem::JacRow jac() const;		// Returns the local Jacobian
                                        // of this equation
    void init(){}
    EqType getType() const {return EqType::Mu_is_m_u;}
    ~Mu(){}
  };
  
  
} // namespace tube




#endif // MU_H
//////////////////////////////////////////////////////////////////////


