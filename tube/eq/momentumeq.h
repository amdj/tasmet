#pragma once

#include "tubeequation.h"
#include "drag.h"

namespace tube{
  SPOILNAMESPACE
  class TubeVertex;
  class Tube;
  class Momentum:public TubeEquation
  {
    const DragResistance* drag=NULL;
  public:
    d Wddt=0,WuL=0,Wu=0,WuR=0;
    d WpL=0,WpR=0;
    Momentum(const TubeVertex& v):TubeEquation(v){}
    virtual void init();
    virtual tasystem::JacRow jac() const;
    virtual enum EqType getType() const { return EqType::Mom;}    
    virtual vd error() const;			// Error in momentum equation at node i
    virtual void show() const;
    virtual void domg(vd& domg) const;

    variable::var momentumFlow() const;
    vd extrapolateMomentumFlow() const;
    tasystem::JacRow dExtrapolateMomentumFlow() const;
  private:
    tasystem::JacCol drhoL() const;
    tasystem::JacCol drho() const;
    tasystem::JacCol drhoR() const;    

    tasystem::JacCol dUL() const; // Are virtual for bc
    tasystem::JacCol dU() const;
    tasystem::JacCol dUR() const;

    tasystem::JacCol dpL() const;
    tasystem::JacCol dpR() const;

  };
}



