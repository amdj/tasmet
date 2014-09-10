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
    
    virtual enum EqType getType() const { return EqType::Mom;}    
    virtual vd error(const TubeVertex&) const;			// Error in momentum equation at node i
    virtual void show() const;
    virtual vd domg(const TubeVertex&) const;
  protected:
    JacCol drhoi(const TubeVertex& v) const;
    JacCol dpi(const TubeVertex&) const;

    JacCol drhoip1(const TubeVertex&) const;
    JacCol dUip1(const TubeVertex&) const;

    JacCol dpim1(const TubeVertex&) const;
    JacCol drhoim1(const TubeVertex&) const; 
    JacCol dpip1(const TubeVertex&) const;
    
    // These two are virtual for the impedance bc's
    virtual JacCol dUim1(const TubeVertex&) const;
    virtual JacCol dUi(const TubeVertex&) const;



    // virtual JacCol dUim2(const TubeVertex&) const;
    // virtual JacCol dUip2(const TubeVertex&) const;
    
    d muR(const TubeVertex&) const;
    d muL(const TubeVertex&) const;
    
  };
}



