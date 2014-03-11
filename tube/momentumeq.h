#pragma once
#include "tubeequation.h"


namespace tube{
  class Momentum:public Equation
  {
  public:
    Momentum(Tube* tube,TubeVertex* gp);
    ~Momentum();
    dmat operator()();	       // Link to Equation operator()
    vd Error();			// Error in momentum equation at node i
    dmat drhoi();
    dmat dUi();

    dmat drhoip1();
    dmat dUip1();

    dmat drhoim1();
    dmat dUim1();

    dmat dpip1();
    dmat dpim1();
    
    
  };
}
