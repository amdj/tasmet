// isentropic.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef ISENTROPIC_H
#define ISENTROPIC_H
#include "tubeequation.h"

namespace tube{
  SPOILNAMESPACE

  class Isentropic:public Equation{
  public:
    Isentropic(const Cell& v):Equation(v){}
    virtual Equation* copy() const {return new Isentropic(*this);}    
    virtual void init();
    virtual enum EqType getType() const { return EqType::Ise;}
    virtual vd error() const;			// Error in Energy equation at node i
    virtual tasystem::JacRow jac() const;
    void show() const;
  private:
    tasystem::JacCol dp() const;
    tasystem::JacCol drho() const;
  };
}

#endif // ISENTROPIC_H
//////////////////////////////////////////////////////////////////////
