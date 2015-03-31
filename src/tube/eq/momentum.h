#pragma once

#include "tubeequation.h"



namespace tube{
  SPOILNAMESPACE
  class Cell;
  class Tube;
  class DragResistance;

  class Momentum:public TubeEquation
  {
    const DragResistance* drag=NULL;
  public:
    d Wddt=0,Wpi=0,Wpim1=0;

    Momentum(const Cell& v):TubeEquation(v){}
    virtual void init();
    virtual tasystem::JacRow jac() const;
    virtual enum EqType getType() const { return EqType::Mom;}    
    virtual vd error() const;			// Error in momentum equation at node i
    virtual void show() const;
    virtual void domg(vd& domg) const;

  private:
    // tasystem::JacCol drho() const;

    // tasystem::JacCol drhoUL() const; // Are virtual for bc
    // tasystem::JacCol drhoULL() const; // Are virtual for bc
    // tasystem::JacCol drhoUR() const;

    // tasystem::JacCol dpi() const;
    // tasystem::JacCol dpim1() const;

  };
}



