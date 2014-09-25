#pragma once

#include "tubeequation.h"
#include "drag.h"

namespace tube{
  SPOILNAMESPACE
  class TubeVertex;
  class Tube;
  class Momentum:public TubeEquation
  {
  public:
    const DragResistance* drag=NULL;
    virtual void init(const Tube& t);
    virtual JacRow jac(const TubeVertex& v) const;
    virtual TubeEquation* copy() const {return new Momentum(*this);}    
    virtual enum EqType getType() const { return EqType::Mom;}    
    virtual vd error(const TubeVertex&) const;			// Error in momentum equation at node i
    virtual void show() const;
    virtual void domg(const TubeVertex&,vd& domg) const;
  protected:
    JacCol drhoim1(const TubeVertex& v) const;
    JacCol drhoi(const TubeVertex& v) const;
    JacCol drhoip1(const TubeVertex& v) const;    

    virtual JacCol dUim1(const TubeVertex&) const; // Are virtual for bc
    virtual JacCol dUi(const TubeVertex&) const;
    virtual JacCol dUip1(const TubeVertex&) const;

    virtual JacCol dpL(const TubeVertex&) const;
    virtual JacCol dpR(const TubeVertex&) const;
    
    // These two are virtual for the impedance bc's





    // virtual JacCol dUim2(const TubeVertex&) const;
    // virtual JacCol dUip2(const TubeVertex&) const;
    
    d muR(const TubeVertex&) const;
    d muL(const TubeVertex&) const;
    
  };
}



