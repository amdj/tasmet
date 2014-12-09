#pragma once
#include "tubeequation.h"

namespace tube{

  class SolidTPrescribed:public TubeEquation
  {
  private:
    const vd* Tsmirror=NULL;
  public:
    SolidTPrescribed(const TubeVertex& v):TubeEquation(v){}
    virtual void init();
    virtual enum EqType getType() const { return EqType::Sol;}
    tasystem::JacRow jac() const;
    
    vd error() const;			// Error in Solidenergy equation at node i
  private:
    tasystem::JacCol dTsi() const;
    
  };
}
