#pragma once
#include "tubeequation.h"

namespace tube{

  class Continuity:public TubeEquation{	// Tube continuity equation 
    d Wddt=0,Wim1=0,Wi=0,Wip1=0;
    // Continuity artificial viscosity weight factor
    d Wart1=0,Wart2=0,Wart3=0,Wart4=0;		       
  public:
    Continuity(const TubeVertex& v):TubeEquation(v){TRACE(15,"Continuity()");}
    virtual void init(const WeightFactors&,const Tube& t);
    virtual enum EqType getType() const { return EqType::Con;}
    virtual void show() const;
    virtual tasystem::JacRow jac() const;
    virtual vd error() const;			// Error in this equation at this node
    virtual void domg(vd&) const;
  private:
    tasystem::JacCol drhoip1() const; // Derivative of continuity equation to density at node
    // i + 1
    tasystem::JacCol drhoi() const;	// Derivative of continuity equation to density at
    // novirtual de i
    tasystem::JacCol dUi() const;	// Derivative of continuity equation to Volume flow at i
    // node (only nonzero for nonconstant grids)
    tasystem::JacCol drhoim1() const;	// Derivative of continuity equation to density at
    // node i - 1
    tasystem::JacCol dUip1() const;	// Derivative of continuity equation to Volume flow
    // at node i + 1
    tasystem::JacCol dUim1() const;	// Derivative of continuity equation to Volume flow
			// at node i - 1
    // virtual JacCol drhoip2() const;
    // virtual JacCol drhoim2() const;

  };				// Continuity class
}				// Namespace tube

