#pragma once

#define MOM_SCALE (1e6)

#include "tubeequation.h"
#include "drag.h"

namespace tube{
  class Momentum:public Equation
  {
  public:
    Momentum(const Tube& tube,TubeVertex& gp);
    ~Momentum();
    d Wuim1,Wui,Wuip1,Wpim1,Wpi,Wpip1; // Weight functions

    vd Error();			// Error in momentum equation at node i
    dmat drhoi();
    dmat dUi();
    dmat dpi();
    dmat drhoip1();
    dmat dUip1();
    dmat dpim1();
    dmat drhoim1();
    dmat dUim1();
    dmat dpip1();

    
    
  };
}
