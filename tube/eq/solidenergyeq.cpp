#include "solidenergyeq.h"
#include "tubevertex.h"

namespace tube{

  SolidTPrescribed::SolidTPrescribed(){
    // Tset is false
  }
  void SolidTPrescribed::setTs(d Tl,d Tr){
    Tset=true;
    this->Tl=Tl;
    this->Tr=Tr;
  }
  void SolidTPrescribed::setTs(d T){
    Tset=true;
    Tl=Tr=T;
  }
  vd SolidTPrescribed::error(const TubeVertex& v) const {		// Error in momentum equation
    TRACE(6,"SolidTPrescribed::Error()");
    // vd error(vertex.gc->Ns,fillwith::zeros);
    assert(v.gc!=NULL);
    assert(v.lg.geom!=NULL);
    vd error;
    if(Tset==false)
      error=v.Ts()-v.gc->T0;
    else
      error=v.Ts()-(Tl+(v.lg.vxi/v.lg.geom->L)*(Tr-Tl));
    return error;
  }
  dmat SolidTPrescribed::dTsi(const TubeVertex& v) const {
    TRACE(0,"SolidTPrescribed:dTsi()");
    // Set solid temperature to zero
    return eye<dmat>(v.gc->Ns,v.gc->Ns);
  }

} // namespace tube
