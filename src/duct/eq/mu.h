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
#include "ductequation.h"

namespace tasystem {
  class JacRow;
} // namespace tasystem

namespace duct {
  SPOILNAMESPACE  
  class Cell;
  class Duct;
  class WeightFactors;

  class MuEq: public Equation{
  public:
    MuEq(const Cell& v):Equation(v){}
    vd error() const;

    tasystem::JacRow jac() const;		// Returns the local Jacobian
                                        // of this equation
    void init(){}
    void show() const;
    EqType getType() const {return EqType::Mu_is_m_u;}
    ~MuEq(){}
  };
  
  
} // namespace duct




#endif // MU_H
//////////////////////////////////////////////////////////////////////


