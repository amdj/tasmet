#include "continuityeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{
  Continuity::Continuity(const Tube& tube,const TubeVertex& gp):
    Equation(tube,gp){
    TRACE(0,"Continuity constructor done");
    wim1=wi=wip1=0;		// Initialize to zero

    // The default boundary implementation is an adiabatic no-slip wall.
    if(i==0){
      wim1=0;
      wi=wRl;
      wip1=WRr;
    }
    else if(i==Ncells-1){
      wi=-wLr;
      wim1=-wLl;
      wip1=0;
    }
    else{ // interior vertex
      wim1=wLl;
      wi=wRl-wLr;
      wip1=wRr;
    }
  }
  vd Continuity::Error(){	// Current error in continuity equation
    vd error(Ns,fillwith::zeros);
    error+=vVf*DDTfd*vertex.rho();
    error+=wi*fDFT*(vertex.rho.tdata()%vertex.U.tdata());
    if(i<Ncells-1){	// Standard implementation of a no-slip (wall)
			// boundary condition
      vd rhoip1=tube.vvertex[i+1].rho.tdata();
      vd Uip1=tube.vvertex[i+1].U.tdata();
      error+=wip1*fDFT*(rhoip1%Uip1);
    }
    if(i>0){ // Standard implementation of a no-slip (wall) boundary condition
      vd rhoim1=tube.vvertex[i-1].rho.tdata();
      vd Uim1=tube.vvertex[i-1].U.tdata();
      error+=wim1*fDFT*(rhoim1%Uim1);
    }
    return error;
  }
  dmat Continuity::drhoi(){
    TRACE(0,"Continuity::drhoi()");
    dmat drhoi=vVf*DDTfd;		// Initialize and add first term
    drhoi+=wi*fDFT*diagtmat(vertex.U)*iDFT;
    return drhoi;
  }
  dmat Continuity::dUi(){
    TRACE(0,"Continuity::dUi()");
    dmat dUi=zero;
    dUi+=wi*fDFT*diagtmat(vertex.rho)*iDFT;
    return dUi;
  }
  dmat Continuity::dUip1(){
    TRACE(0,"Continuity::dUip1()");
    dmat dUip1=zero;
    if(i<Ncells-1)
      dUip1+=wip1*fDFT*diagtmat(tube.vvertex[i+1].rho)*iDFT;
    return dUip1;
  }
  dmat Continuity::dUim1(){
    TRACE(0,"Continuity::dUim1()");
    dmat dUim1=zero;
    if(i>0)
      dUim1+=wim1*fDFT*diagtmat(tube.vvertex[i-1].rho)*iDFT;
    return dUim1;
  }
  dmat Continuity::drhoip1(){
    TRACE(0,"Continuity::drhoip1()");
    dmat drhoip1=zero;
    if(i<Ncells-1)
      drhoip1=wip1*fDFT*diagtmat(tube.vvertex[i+1].U)*iDFT;
    return drhoip1;
  }
  dmat Continuity::drhoim1(){
    TRACE(0,"Continuity::drhoim1()");
    dmat drhoim1=zero;
    if(i>0)
      drhoim1+=wim1*fDFT*diagtmat(tube.vvertex[i-1].U)*iDFT;
    return drhoim1;
  }

  Continuity::~Continuity(){}

} // Namespace tube
