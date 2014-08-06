#pragma once
#ifndef _ISENTROPICEQ_H_
#define _ISENTROPICEQ_H_




#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE

  class Isentropic:public TubeEquation{
  public:
    virtual enum EqType getType() const { return EqType::Ise;}
    virtual vd error(const TubeVertex&) const;			// Error in Energy equation at node i
    virtual dmat dpi(const TubeVertex&) const;
    virtual dmat drhoi(const TubeVertex&) const;
  };
}

#endif /* _ISENTROPICEQ_H_ */
