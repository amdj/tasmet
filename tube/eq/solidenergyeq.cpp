#include "solidenergyeq.h"
#include "tubevertex.h"

namespace tube{

  SolidTPrescribed::SolidTPrescribed(){
    // Tset is false
  }
  void SolidTPrescribed::init(const Tube& t) {
    TRACE(6,"SolidTPrescribed::init(t)");
  }
  void SolidTPrescribed::setTs(const Geom& g,d Tl,d Tr){
    Tset=true;
    this->vx=g.vx;
    this->L=g.L;
    this->Tl=Tl;
    this->Tr=Tr;
  }
  void SolidTPrescribed::setTs(const Geom& g,d T){
    setTs(g,T,T);
  }
  vd SolidTPrescribed::error(const TubeVertex& v) const {		// Error in momentum equation
    TRACE(6,"SolidTPrescribed::Error()");
    // vd error(vertex.gc->Ns,fillwith::zeros);
    assert(v.gc!=NULL);
    vd error;
    if(Tset==false){
      error=v.Ts();
      error(0)+=-v.gc->T0;
    }
    else{
      error=v.Ts();
      error(0)+=-(Tl+(vx(v.i)/L)*(Tr-Tl));
    }
    return error;
  }
  dmat SolidTPrescribed::dTsi(const TubeVertex& v) const {
    TRACE(0,"SolidTPrescribed:dTsi()");
    // Set solid temperature to zero
    return eye<dmat>(v.gc->Ns,v.gc->Ns);
  }

} // namespace tube
