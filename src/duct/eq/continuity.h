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
#include "ductequation.h"

namespace duct{

  class Continuity:public Equation{	// Duct continuity equation 
    d Wddt=0;
  public:
    Continuity(const Cell& v):Equation(v){TRACE(15,"Continuity()");}
    ~Continuity();
    virtual void init();
    virtual enum EqType getType() const { return EqType::Con;}
    virtual void show() const;
    virtual tasystem::JacRow jac() const;
    virtual vd error() const;			// Error in this equation at this node
    virtual void domg(vd&) const;
    // vd massFlow() const;
    static vd extrapolateMassFlow(const Cell&);
    static tasystem::JacRow dExtrapolateMassFlow(const Cell&);
  };				// Continuity class
  
} // Namespace duct

#endif // CONTINUITY_H
//////////////////////////////////////////////////////////////////////
