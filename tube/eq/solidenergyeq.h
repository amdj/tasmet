#pragma once
#include "tubeequation.h"


namespace tube{

  class SolidTPrescribed:public TubeEquation
  {
  private:
    bool Tset=false;
    d Tl=0,Tr=0;
  public:
    SolidTPrescribed();
    void setTs(d Tl,d Tr);
    void setTs(d T);
    virtual enum EqType getType() const { return EqType::Sol;}
    vd error(const TubeVertex&) const;			// Error in Solidenergy equation at node i
    dmat dTsi(const TubeVertex&) const;
    
  };
}
