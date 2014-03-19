#pragma once
#include "tubeequation.h"


namespace tube{
  class Energy:public Equation
  {
  public:
    Energy(const Tube& tube,const TubeVertex& gp);
    ~Energy();
    vd Error();			// Error in Energy equation at node i
    dmat dpi();
    dmat dUi();
    dmat dTi();
    
    dmat dpip1();
    dmat dUip1();
    dmat dTip1();
    
    dmat dUim1();
    dmat dpim1();
    dmat dTim1();
    
    vd kappaL();		// Thermal conducticity at the left
				// boundary of the vertex
    vd kappaR();		// Thermal conductivity at the right
				// boundary of the vertex
    d gamma();			// Time-avg ratio of specific heats
  };
}
