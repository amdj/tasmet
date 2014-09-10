#pragma once
#ifndef _ISENTROPICEQ_H_
#define _ISENTROPICEQ_H_




#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE

  class Isentropic:public TubeEquation{
  public:
    virtual void init(const Tube& t);    
    virtual enum EqType getType() const { return EqType::Ise;}
    virtual vd error(const TubeVertex&) const;			// Error in Energy equation at node i
    virtual JacRow jac(const TubeVertex&) const;
  private:
    JacCol dpi(const TubeVertex&) const;
    JacCol drhoi(const TubeVertex&) const;
  };
}

#endif /* _ISENTROPICEQ_H_ */
