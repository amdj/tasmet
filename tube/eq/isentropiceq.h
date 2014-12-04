#pragma once
#ifndef _ISENTROPICEQ_H_
#define _ISENTROPICEQ_H_




#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE

  class Isentropic:public TubeEquation{
  public:
    Isentropic(const TubeVertex& v):TubeEquation(v){}
    virtual TubeEquation* copy() const {return new Isentropic(*this);}    
    virtual void init(const WeightFactors&,const Tube&);    
    virtual enum EqType getType() const { return EqType::Ise;}
    virtual vd error(const TubeVertex&) const;			// Error in Energy equation at node i
    virtual tasystem::JacRow jac(const TubeVertex&) const;
  private:
    virtual tasystem::JacCol dpL(const TubeVertex&) const;
    virtual tasystem::JacCol dpR(const TubeVertex&) const;    
    tasystem::JacCol drhoi(const TubeVertex&) const;
  };
}

#endif /* _ISENTROPICEQ_H_ */
