#pragma once
#include "tubeequation.h"
#include "geom.h"

namespace tube{
  using segment::Geom;
  class SolidTPrescribed:public TubeEquation
  {
  private:
    vd Tsmirror;
  public:
    void  setTMirror(vd Ts) {Tsmirror=Ts;}
    SolidTPrescribed();
    void init(const Tube& t);
    virtual enum EqType getType() const { return EqType::Sol;}
    vd error(const TubeVertex&) const;			// Error in Solidenergy equation at node i
    dmat dTsi(const TubeVertex&) const;
    
  };
}
