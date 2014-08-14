#pragma once
#include "tubeequation.h"
#include "geom.h"

namespace tube{
  using segment::Geom;
  class SolidTPrescribed:public TubeEquation
  {
  private:
    bool Tset=false;
    d Tl=0,Tr=0;
    vd xv;
    d L;
  public:
    SolidTPrescribed();
    void init(const Tube& t);
    void setTs(const Geom&,d Tl,d Tr);
    void setTs(const Geom&,d T);
    virtual enum EqType getType() const { return EqType::Sol;}
    vd error(const TubeVertex&) const;			// Error in Solidenergy equation at node i
    dmat dTsi(const TubeVertex&) const;
    
  };
}
