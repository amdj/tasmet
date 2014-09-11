#pragma once
#include "tubeequation.h"

namespace tube{



  class Continuity:public TubeEquation{	// Tube continuity equation 
  public:
    virtual enum EqType getType() const { return EqType::Con;}
    virtual void show() const;
    virtual void init(const Tube& t);
    virtual JacRow jac(const TubeVertex&) const;
    virtual vd error(const TubeVertex&) const;			// Error in this equation at this node
    virtual vd domg(const TubeVertex&) const;
    virtual TubeEquation* copy() const {return new Continuity(*this);}
  private:
    JacCol drhoip1(const TubeVertex&) const; // Derivative of continuity equation to density at node
    // i + 1
    JacCol drhoi(const TubeVertex&) const;	// Derivative of continuity equation to density at
    // novirtual de i
    JacCol dUi(const TubeVertex&) const;	// Derivative of continuity equation to Volume flow at i
    // node (const TubeVertex&only nonzero for nonconstant grids)
    JacCol drhoim1(const TubeVertex&) const;	// Derivative of continuity equation to density at
    // node i - 1
     JacCol dUip1(const TubeVertex&) const;	// Derivative of continuity equation to Volume flow
    // at node i + 1
    JacCol dUim1(const TubeVertex&) const;	// Derivative of continuity equation to Volume flow
			// at node i - 1
    // virtual JacCol drhoip2(const TubeVertex&) const;
    // virtual JacCol drhoim2(const TubeVertex&) const;

  };				// Continuity class
}				// Namespace tube

