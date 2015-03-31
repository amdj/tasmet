#pragma once
#ifndef _ISENTROPICEQ_H_
#define _ISENTROPICEQ_H_




#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE

  class Isentropic:public TubeEquation{
  public:
    Isentropic(const Cell& v):TubeEquation(v){}
    virtual TubeEquation* copy() const {return new Isentropic(*this);}    
    virtual void init();
    virtual enum EqType getType() const { return EqType::Ise;}
    virtual vd error() const;			// Error in Energy equation at node i
    virtual tasystem::JacRow jac() const;
  private:
    tasystem::JacCol dp() const;
    tasystem::JacCol drho() const;
  };
}

#endif /* _ISENTROPICEQ_H_ */
