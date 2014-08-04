#pragma once
#ifndef _ISENTROPICEQ_H_
#define _ISENTROPICEQ_H_




#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE

  class Isentropic:public TubeEquation{
  public:
    virtual enum EqType getType() const { return EqType::Ise;}
    Isentropic();
    ~Isentropic();
    vd error(const TubeVertex&) const;			// Error in Energy equation at node i
    dmat dpi(const TubeVertex&) const;
    dmat drhoi(const TubeVertex&) const;

  };
}

#endif /* _ISENTROPICEQ_H_ */
