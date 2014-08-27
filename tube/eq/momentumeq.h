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
    virtual vd domg(const TubeVertex&) const;
    virtual dmat drhoi(const TubeVertex& v) const;
    virtual dmat dUi(const TubeVertex&) const;
    virtual dmat dpi(const TubeVertex&) const;
    virtual dmat drhoip1(const TubeVertex&) const;
    virtual dmat dUip1(const TubeVertex&) const;
    virtual dmat dpim1(const TubeVertex&) const;
    virtual dmat drhoim1(const TubeVertex&) const;
    virtual dmat dUim1(const TubeVertex&) const;
    virtual dmat dpip1(const TubeVertex&) const;

    virtual dmat dUim2(const TubeVertex&) const;
    virtual dmat dUip2(const TubeVertex&) const;
    
    d muR(const TubeVertex&) const;
    d muL(const TubeVertex&) const;
    
  };
}



