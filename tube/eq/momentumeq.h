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
    virtual enum EqType getType() const { return EqType::Mom;}    
    virtual vd error(const TubeVertex&) const;			// Error in momentum equation at node i
    virtual void show() const;
    dmat drhoi(const TubeVertex& v) const;
    dmat dUi(const TubeVertex&) const;
    dmat dpi(const TubeVertex&) const;
    dmat drhoip1(const TubeVertex&) const;
    dmat dUip1(const TubeVertex&) const;
    dmat dpim1(const TubeVertex&) const;
    dmat drhoim1(const TubeVertex&) const;
    dmat dUim1(const TubeVertex&) const;
    dmat dpip1(const TubeVertex&) const;

    dmat dUim2(const TubeVertex&) const;
    dmat dUip2(const TubeVertex&) const;
    
    d muR(const TubeVertex&) const;
    d muL(const TubeVertex&) const;
    
  };
}



