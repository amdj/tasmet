#pragma once
#include "tubeequation.h"

namespace tube{

  class SolidTPrescribed:public TubeEquation
  {
  private:
    const vd* Tsmirror=NULL;
  public:
    SolidTPrescribed();
    void init(const Tube& t);
    virtual enum EqType getType() const { return EqType::Sol;}
    JacRow jac(const TubeVertex&) const;
    
    vd error(const TubeVertex&) const;			// Error in Solidenergy equation at node i
  private:
    JacCol dTsi(const TubeVertex&) const;
    
  };
}
