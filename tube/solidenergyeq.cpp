#include "solidenergyeq.h"
#include "tubevertex.h"

namespace tube{

  Solidenergy::Solidenergy(TubeVertex& gp):TubeEquation(gp){
    TRACE(7,"SolidSolidenergy constructor done");
  }

  vd Solidenergy::Error(){		// Error in momentum equation
    TRACE(6,"Solidenergy::Error()");
    // vd error(gc->Ns,fillwith::zeros);
    vd error=vertex.Ts();
    return error;
  }
  dmat Solidenergy::dTsi(){
    TRACE(0,"Solidenergy:dTsi()");
    // Set solid temperature to zero
    return eye<dmat>(gc->Ns,gc->Ns);
  }
  Solidenergy::~Solidenergy(){
    TRACE(-5,"Solidenergy destructor");
  }
} // namespace tube
