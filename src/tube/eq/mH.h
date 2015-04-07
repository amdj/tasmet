// mHEq.h
//
// Author: J.A. de Jong 
//
// Description:
// This equation relates the total enthalpy flow at the left boundary
// to the variables at the vertex. This equation should only be used
// for interior nodes.
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef MH_H
#define MH_H
#include "vtypes.h"
#include "tubeequation.h"

namespace tasystem {
  class JacRow;
} // namespace tasystem

namespace tube {
  SPOILNAMESPACE  
  class Cell;
  class Tube;

  class mHEq: public Equation{
    d WLr=0,WLl=0;
  public:
    mHEq(const Cell& v):Equation(v){}
    vd error() const;

    tasystem::JacRow jac() const;		// Returns the local Jacobian
                                        // of this equation
    void init();
    void show() const;
    EqType getType() const {return EqType::mH_is_m_H;}
    ~mHEq(){}
  };
  
  
} // namespace tube



#endif // MH_H
//////////////////////////////////////////////////////////////////////


