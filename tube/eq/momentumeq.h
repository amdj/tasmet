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
    d Wart1=0,Wart2=0,Wart3=0,Wart4=0;		       
    d Wddt=0,Wuim1=0,Wui=0,Wuip1=0;
    d WpL=0,WpR=0;
    Momentum(const TubeVertex& v):TubeEquation(v){}
    virtual void init(const WeightFactors&,const Tube& t);
    virtual tasystem::JacRow jac() const;
    virtual enum EqType getType() const { return EqType::Mom;}    
    virtual vd error() const;			// Error in momentum equation at node i
    virtual void show() const;
    virtual void domg(vd& domg) const;
  protected:
    tasystem::JacCol drhoim1() const;
    tasystem::JacCol drhoi() const;
    tasystem::JacCol drhoip1() const;    

    virtual tasystem::JacCol dUim1() const; // Are virtual for bc
    virtual tasystem::JacCol dUi() const;
    virtual tasystem::JacCol dUip1() const;

    virtual tasystem::JacCol dpL() const;
    virtual tasystem::JacCol dpR() const;
    
    // These two are virtual for the impedance bc's





    // virtual JacCol dUim2() const;
    // virtual JacCol dUip2() const;
    
    d muR() const;
    d muL() const;
    
  };
}



