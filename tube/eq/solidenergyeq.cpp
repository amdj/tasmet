#include "solidenergyeq.h"
#include "tubevertex.h"

namespace tube{

  vd Solidenergy::error(const TubeVertex& v) const {		// Error in momentum equation
    TRACE(6,"Solidenergy::Error()");
    // vd error(vertex.gc->Ns,fillwith::zeros);
    vd error=v.Ts();
    return error;
  }
  dmat Solidenergy::dTsi(const TubeVertex& v) const {
    TRACE(0,"Solidenergy:dTsi()");
    // Set solid temperature to zero
    return eye<dmat>(v.gc->Ns,v.gc->Ns);
  }

} // namespace tube
