#pragma once
#include "tubeequation.h"

namespace tube{

  class SolidTPrescribed:public TubeEquation
  {
  private:
    const vd* Tsmirror=NULL;
  public:
    SolidTPrescribed(const TubeVertex& v):TubeEquation(v){}
    void init(const WeightFactors&,const Tube&);
    virtual enum EqType getType() const { return EqType::Sol;}
    tasystem::JacRow jac(const TubeVertex&) const;
    
    vd error(const TubeVertex&) const;			// Error in Solidenergy equation at node i
  private:
    tasystem::JacCol dTsi(const TubeVertex&) const;
    
  };
}
