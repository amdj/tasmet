#pragma once
#include "tubeequation.h"


namespace tube{
  class Energy:public Equation
  {
  public:
    Energy(const Tube& tube,const TubeVertex& gp);
    ~Energy();
    dmat operator()();
    vd Error();			// Error in Energy equation at node i
    dmat dpi();
    dmat dUi();

    dmat dUip1();
    dmat dUim1();

    dmat dpip1();
    dmat dpim1();
    
    
  };
}
