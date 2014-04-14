#include "continuityeq.h"
#include "vertex.h"
#include "tube.h"

#define CONT_VISCOSITY
#define CONT_SCALE (1.0)//(pow(gc.c0,2)) //(pow(gc.c0,2)/gc.omg)
#define CONT_SCALE0 (1.0)//(1.0/pow(gc.M,2))



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
    TRACE(5,"Continuity zero scale:"<<CONT_SCALE0);
  }
  vd Continuity::Error(){	// Current error in continuity equation
    vd error(Ns,fillwith::zeros);
    error+=vVf*DDTfd*vertex.rho();
    error+=Wi*fDFT*(vertex.rho.tdata()%vertex.U.tdata());
    
    if(i>0){
      // Standard implementation of a no-slip (wall) boundary
     // condition
      vd rhoim1=left->rho.tdata();
      vd Uim1=left->U.tdata();
      error+=Wim1*fDFT*(rhoim1%Uim1);
    }
    if(i<Ncells-1){
      // Standard implementation of a no-slip (wall) boundary
      // condition
      vd rhoip1=right->rho.tdata();
      vd Uip1=right->U.tdata();
      error+=Wip1*fDFT*(rhoip1%Uip1);
    }

    // Artificial viscosity term
    // TRACE(10,D_l());
    // TRACE(10,D_r());
    
    #ifdef CONT_VISCOSITY
    if(i>0 && i<Ncells-1){
      error+=-D_r()*(right->rho() -vertex.rho())*vSf;
      error+= D_l()*(vertex.rho() -left->rho())  *vSf;
    }
    else if(i==0){		// First vertex
      error+=-D_r()*(right->right->rho()-right->rho())*vSf;
      error+=D_l()*(right->rho()-vertex.rho())*vSf;
    }
    else {			// Last vertex
      error+=-D_r()*(vertex.rho()-left->rho())*vSf;
      error+=D_l()*(left->rho()-left->left->rho())*vSf;
    }
    #endif
    
    
    // (Boundary) source term
    error+=vertex.csource();
    error(0)*=CONT_SCALE0;
    return CONT_SCALE*error;
  }
  dmat Continuity::drhoi(){
    TRACE(0,"Continuity::drhoi()");
    dmat drhoi=vVf*DDTfd;		// Initialize and add first term
    drhoi+=Wi*fDFT*diagtmat(vertex.U)*iDFT;

    // Artificial viscosity terms
    #ifdef CONT_VISCOSITY
    if(i>0 && i<Ncells-1){
      drhoi+=(D_l()+D_r())*vSf;	// Middle vertex
    }
    else if(i==0)
      drhoi+=-D_l()*vSf;	// First vertex
    else		
      drhoi+=-D_r()*vSf;	// Last vertex
    #endif
    drhoi.row(0)*=CONT_SCALE0;    
    return CONT_SCALE*drhoi;
  }
  dmat Continuity::dUi(){
    TRACE(0,"Continuity::dUi()");
    dmat dUi=zero;
    dUi+=Wi*fDFT*diagtmat(vertex.rho)*iDFT;
    dUi.row(0)*=CONT_SCALE0;    
    return CONT_SCALE*dUi;
  }
  dmat Continuity::dUip1(){
    TRACE(0,"Continuity::dUip1()");
    dmat dUip1=zero;
    if(i<Ncells-1)
      dUip1+=Wip1*fDFT*diagtmat(right->rho)*iDFT;
    dUip1.row(0)*=CONT_SCALE0;
    return CONT_SCALE*dUip1;
  }
  dmat Continuity::dUim1(){
    TRACE(0,"Continuity::dUim1()");
    dmat dUim1=zero;
    if(i>0)
      dUim1+=Wim1*fDFT*diagtmat(left->rho)*iDFT;
    dUim1.row(0)*=CONT_SCALE0;
    return CONT_SCALE*dUim1;
  }
  dmat Continuity::drhoip1(){
    TRACE(0,"Continuity::drhoip1()");
    dmat drhoip1=zero;
    if(i<Ncells-1)
      drhoip1=Wip1*fDFT*diagtmat(right->U)*iDFT;

    // Artificial viscosity terms
    #ifdef CONT_VISCOSITY    
    if(i>0 && i<Ncells-1){
      drhoip1+=-D_r()*vSf;
    }
    else if(i==0){		// First vertex
      drhoip1+=(D_l()+D_r())*vSf;
    }
    #endif
    
    drhoip1.row(0)*=CONT_SCALE0;
    return CONT_SCALE*drhoip1;
  }
  dmat Continuity::drhoim1(){
    TRACE(0,"Continuity::drhoim1()");
    dmat drhoim1=zero;
    if(i>0)
      drhoim1+=Wim1*fDFT*diagtmat(left->U)*iDFT;

    // Artificial viscosity terms
    #ifdef CONT_VISCOSITY
    if(i>0 && i<Ncells-1){
      drhoim1+=-D_l()*vSf;
    }
    else if(i==Ncells-1){		// Last vertex
      drhoim1+=(D_l()+D_r())*vSf;
    }
    #endif

    drhoim1.row(0)*=CONT_SCALE0;
    return CONT_SCALE*drhoim1;
  }
  dmat Continuity::drhoim2(){
    dmat drhoim2=zero;
    if(i==Ncells-1)
      drhoim2+=-D_l()*vSf;
    return drhoim2;
  }
  dmat Continuity::drhoip2(){
    dmat drhoip2=zero;
    if(i==0)
      drhoip2+=-D_r()*vSf;
    return drhoip2;

  }
  Continuity::~Continuity(){}

} // Namespace tube

