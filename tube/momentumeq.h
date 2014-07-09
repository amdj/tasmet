#pragma once

#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE
  class TubeVertex;
  
  class Momentum:public TubeEquation
  {
  public:
    Momentum(TubeVertex& gp);
    ~Momentum();
    d Wuim1,Wui,Wuip1,Wpim1,Wpi,Wpip1; // Weight functions
    d Wddt;			       // Time-derivative weight function
    virtual vd Error();			// Error in momentum equation at node i
    dmat drhoi();
    dmat dUi();
    dmat dpi();
    dmat drhoip1();
    dmat dUip1();
    dmat dpim1();
    dmat drhoim1();
    dmat dUim1();
    dmat dpip1();

    dmat dUim2();
    dmat dUip2();
    
    d muR();
    d muL();
    
  };
}



