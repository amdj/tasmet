#pragma once
#include "tubeequation.h"

namespace tube{

  class Continuity:public TubeEquation{	// Tube continuity equation 
    d Wddt=0,WL=0,Wi=0,WR=0;
  public:
    Continuity(const Cell& v):TubeEquation(v){TRACE(15,"Continuity()");}
    virtual void init();
    virtual enum EqType getType() const { return EqType::Con;}
    virtual void show() const;
    virtual tasystem::JacRow jac() const;
    virtual vd error() const;			// Error in this equation at this node
    virtual void domg(vd&) const;
    // vd massFlow() const;

  private:
    tasystem::JacCol drho() const;	// Derivative of continuity equation to density at
    // novirtual de i
    tasystem::JacCol drhoUR() const;	// Derivative of continuity equation to Volume flow at i
    // node (only nonzero for nonconstant grids)
    tasystem::JacCol drhoUL() const;	// Derivative of continuity equation to density at
  };				// Continuity class
}				// Namespace tube

