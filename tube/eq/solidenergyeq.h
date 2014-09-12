#pragma once
#include "tubeequation.h"
#include "geom.h"

namespace tube{
  using segment::Geom;
  class SolidTPrescribed:public TubeEquation
  {
  private:
    const vd* Tsmirror=NULL;
  public:
    virtual TubeEquation* copy() const {return new SolidTPrescribed(*this);}    
    SolidTPrescribed();
    void init(const Tube& t);
    virtual enum EqType getType() const { return EqType::Sol;}
    JacRow jac(const TubeVertex&) const;
    
    vd error(const TubeVertex&) const;			// Error in Solidenergy equation at node i
  private:
    JacCol dTsi(const TubeVertex&) const;
    
  };
}
