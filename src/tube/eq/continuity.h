// continuity.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef CONTINUITY_H
#define CONTINUITY_H
#include "tubeequation.h"

namespace tube{

  class Continuity:public Equation{	// Tube continuity equation 
    d Wddt=0;
  public:
    Continuity(const Cell& v):Equation(v){TRACE(15,"Continuity()");}
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
    tasystem::JacCol dmR() const;	// Derivative of continuity equation to Volume flow at i
    // node (only nonzero for nonconstant grids)
    tasystem::JacCol dmL() const;	// Derivative of continuity equation to density at
  };				// Continuity class
}				// Namespace tube

#endif // CONTINUITY_H
//////////////////////////////////////////////////////////////////////
