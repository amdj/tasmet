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
    ~Continuity();
    virtual void init();
    virtual enum EqType getType() const { return EqType::Con;}
    virtual void show() const;
    virtual tasystem::JacRow jac() const;
    virtual vd error() const;			// Error in this equation at this node
    virtual void domg(vd&) const;
    // vd massFlow() const;
    static vd extrapolateMassFlow(const Cell&);
    static vd extrapolateDensity(const Cell&);
    static tasystem::JacRow dExtrapolateMassFlow(const Cell&);
    static tasystem::JacRow dExtrapolateDensity(const Cell&);
  };				// Continuity class
  
} // Namespace tube

#endif // CONTINUITY_H
//////////////////////////////////////////////////////////////////////
