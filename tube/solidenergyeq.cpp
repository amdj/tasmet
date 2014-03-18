#include "solidenergyeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{

  Solidenergy::Solidenergy(const Tube& tube,const TubeVertex& gp):Equation(tube,gp){
    TRACE(0,"SolidSolidenergy constructor done");
  }
  dmat Solidenergy::operator()(){
    TRACE(0,"Solidenergy::operator()");
    return Equation::operator()();
  }
  vd Solidenergy::Error(){		// Error in momentum equation
    TRACE(0,"Solidenergy::Error()");
    vd error(Ns,fillwith::zeros);
    error=vertex->Ts();
    return error;
  }
  dmat Solidenergy::dTsi(){
    TRACE(0,"Solidenergy:dTsi()");
    // Set solid temperature to zero
    return eye<dmat>(Ns,Ns);
  }
  Solidenergy::~Solidenergy(){
    TRACE(-5,"Solidenergy destructor");
  }
} // namespace tube
