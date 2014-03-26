#include "continuityeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{
  Continuity::Continuity(const Tube& tube,TubeVertex& gp):
    Equation(tube,gp){
    TRACE(0,"Continuity constructor done");
    Wim1=Wi=Wip1=0;		// Initialize to zero

    // The default boundary implementation is an adiabatic no-slip wall.
    if(i==0){
      Wim1=0;
      Wi=wRl;
      Wip1=wRr;
    }
    else if(i==Ncells-1){
      Wi=-wLr;
      Wim1=-wLl;
      Wip1=0;
    }
    else{ // interior vertex
      Wim1=-wLl;
      Wi=wRl-wLr;
      Wip1=wRr;
    }
  }
  vd Continuity::Error(){	// Current error in continuity equation
    vd error(Ns,fillwith::zeros);
    error+=vVf*DDTfd*vertex.rho();
    error+=Wi*fDFT*(vertex.rho.tdata()%vertex.U.tdata());
    if(i<Ncells-1){	// Standard implementation of a no-slip (wall)
			// boundary condition
      vd rhoip1=right->rho.tdata();
      vd Uip1=right->U.tdata();
      error+=Wip1*fDFT*(rhoip1%Uip1);
    }
    if(i>0){ // Standard implementation of a no-slip (wall) boundary condition
      vd rhoim1=left->rho.tdata();
      vd Uim1=left->U.tdata();
      error+=Wim1*fDFT*(rhoim1%Uim1);
    }
    // (Boundary) source term
    error+=vertex.csource();
    return error;
  }
  dmat Continuity::drhoi(){
    TRACE(0,"Continuity::drhoi()");
    dmat drhoi=vVf*DDTfd;		// Initialize and add first term
    drhoi+=Wi*fDFT*diagtmat(vertex.U)*iDFT;
    return drhoi;
  }
  dmat Continuity::dUi(){
    TRACE(0,"Continuity::dUi()");
    dmat dUi=zero;
    dUi+=Wi*fDFT*diagtmat(vertex.rho)*iDFT;
    return dUi;
  }
  dmat Continuity::dUip1(){
    TRACE(0,"Continuity::dUip1()");
    dmat dUip1=zero;
    if(i<Ncells-1)
      dUip1+=Wip1*fDFT*diagtmat(right->rho)*iDFT;
    return dUip1;
  }
  dmat Continuity::dUim1(){
    TRACE(0,"Continuity::dUim1()");
    dmat dUim1=zero;
    if(i>0)
      dUim1+=Wim1*fDFT*diagtmat(left->rho)*iDFT;
    return dUim1;
  }
  dmat Continuity::drhoip1(){
    TRACE(0,"Continuity::drhoip1()");
    dmat drhoip1=zero;
    if(i<Ncells-1)
      drhoip1=Wip1*fDFT*diagtmat(right->U)*iDFT;
    return drhoip1;
  }
  dmat Continuity::drhoim1(){
    TRACE(0,"Continuity::drhoim1()");
    dmat drhoim1=zero;
    if(i>0)
      drhoim1+=Wim1*fDFT*diagtmat(left->U)*iDFT;
    return drhoim1;
  }

  Continuity::~Continuity(){}

} // Namespace tube
