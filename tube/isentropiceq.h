#pragma once
#ifndef _ISENTROPICEQ_H_
#define _ISENTROPICEQ_H_




#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE

  class Isentropic:public TubeEquation{
  public:
    Isentropic(TubeVertex& gp);
    ~Isentropic();
    vd Error();			// Error in Energy equation at node i
    dmat dpi();
    dmat drhoi();

  };
}

#endif /* _ISENTROPICEQ_H_ */
