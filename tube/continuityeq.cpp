#include "continuityeq.h"
#include "tubevertex.h"

#define CONT_VISCOSITY
#define CONT_SCALE (1.0)//(pow(v.gc.c0,2)) //(pow(v.gc.c0,2)/v.gc.omg)
#define CONT_SCALE0 (1.0)//(1.0/pow(v.gc.M,2))



namespace tube{
  Continuity::Continuity(TubeVertex& gp):
    TubeEquation(gp){
    TRACE(6,"Continuity constructor done");
    Wim1=Wi=Wip1=Wddt=0;		// Initialize to zero
  }

  vd Continuity::Error(){	// Current error in continuity equation
    // The default boundary implementation is an adiabatic no-slip wall.
    TRACE(6,"Continuity::Error()");
    vd error(v.gc->Ns,fillwith::zeros);
    error+=Wddt*v.gc->DDTfd*v.rho();
    error+=Wi*v.gc->fDFT*(v.rho.tdata()%v.U.tdata());
    if(v.i>0 || (v.i==0 && v.left!=NULL)){
      // Standard implementation of a no-slip (wall) boundary
      // condition
      vd rhoim1=v.left->rho.tdata();
      vd Uim1=v.left->U.tdata();
      error+=Wim1*v.gc->fDFT*(rhoim1%Uim1);
    }
    if(v.i<v.Ncells-1 || (v.i==v.Ncells-1 && v.right!=NULL) ){
      // Standard implementation of a no-slip (wall) boundary
      // condition
      vd rhoip1=v.right->rho.tdata();
      vd Uip1=v.right->U.tdata();
      error+=Wip1*v.gc->fDFT*(rhoip1%Uip1);
    }

    #ifdef CONT_VISCOSITY
    const d& vVf=v.vVf;
    if(v.i>0 && v.i<v.Ncells-1){
      error+=-d_r()*(v.right->rho() -v.rho())*vVf;
      error+= d_l()*(v.rho() -v.left->rho())  *vVf;
    }
    else if(v.i==0){		// First v
      error+=-d_r()*(v.right->right->rho()-v.right->rho())*vVf;
      error+=d_l()*(v.right->rho()-v.rho())*vVf;
    }
    else {			// Last v
      error+=-d_r()*(v.rho()-v.left->rho())*vVf;
      error+=d_l()*(v.left->rho()-v.left->left->rho())*vVf;
    }
    #endif
    
    TRACE(6,"Continuity::Error()");    
    // (Boundary) source term
    error+=v.csource();
    TRACE(6,"Continuity::Error()");
    error(0)*=CONT_SCALE0;
    return CONT_SCALE*error;
  }
  dmat Continuity::drhoi(){
    TRACE(0,"Continuity::drhoi()");
    dmat drhoi=Wddt*v.gc->DDTfd;		// Initialize and add first term
    drhoi+=Wi*v.gc->fDFT*diagtmat(v.U)*v.gc->iDFT;

    // Artificial viscosity terms
    #ifdef CONT_VISCOSITY
    const d& vVf=v.vVf;
    const d& vSf=v.vSf;
    if(v.i>0 && v.i<v.Ncells-1){
      drhoi+=(d_l()+d_r())*vVf;	// Middle vertex
    }
    else if(v.i==0)
      drhoi+=-d_l()*vVf;	// First vertex
    else		
      drhoi+=-d_r()*vVf;	// Last vertex
    #endif
    drhoi.row(0)*=CONT_SCALE0;    
    return CONT_SCALE*drhoi;
  }
  dmat Continuity::dUi(){
    TRACE(0,"Continuity::dUi()");
    dmat dUi=zero;
    dUi+=Wi*v.gc->fDFT*diagtmat(v.rho)*v.gc->iDFT;
    dUi.row(0)*=CONT_SCALE0;    
    return CONT_SCALE*dUi;
  }
  dmat Continuity::dUip1(){
    TRACE(0,"Continuity::dUip1()");
    dmat dUip1=zero;

    if(v.right!=NULL)
      dUip1+=Wip1*v.gc->fDFT*diagtmat(v.right->rho)*v.gc->iDFT;
    dUip1.row(0)*=CONT_SCALE0;
    return CONT_SCALE*dUip1;
  }
  dmat Continuity::dUim1(){
    TRACE(0,"Continuity::dUim1()");

    dmat dUim1=zero;
    if(v.left!=NULL)
      dUim1+=Wim1*v.gc->fDFT*diagtmat(v.left->rho)*v.gc->iDFT;
    dUim1.row(0)*=CONT_SCALE0;
    return CONT_SCALE*dUim1;
  }
  dmat Continuity::drhoip1(){
    TRACE(0,"Continuity::drhoip1()");
    dmat drhoip1=zero;

    if(v.i<v.Ncells-1 || v.right!=NULL)
      drhoip1=Wip1*v.gc->fDFT*diagtmat(v.right->U)*v.gc->iDFT;

    // Artificial viscosity terms
    #ifdef CONT_VISCOSITY    
    const d& vVf=v.vVf;
    if(v.i>0 && v.i<v.Ncells-1){
      drhoip1+=-d_r()*vVf;
    }
    else if(v.i==0){		// First vertex
      drhoip1+=(d_l()+d_r())*vVf;
    }
    #endif
    
    drhoip1.row(0)*=CONT_SCALE0;
    return CONT_SCALE*drhoip1;
  }
  dmat Continuity::drhoim1(){
    TRACE(0,"Continuity::drhoim1()");

    dmat drhoim1=zero;
    if(v.left!=NULL)
      drhoim1+=Wim1*v.gc->fDFT*diagtmat(v.left->U)*v.gc->iDFT;

    // Artificial viscosity terms
    #ifdef CONT_VISCOSITY
    const d& vVf=v.vVf;
    if(v.i>0 && v.i<v.Ncells-1){
      drhoim1+=-d_l()*vVf;
    }
    else if(v.i==v.Ncells-1){		// Last vertex
      drhoim1+=(d_l()+d_r())*vVf;
    }
    #endif

    drhoim1.row(0)*=CONT_SCALE0;
    return CONT_SCALE*drhoim1;
  }
  dmat Continuity::drhoim2(){
    dmat drhoim2=zero;
    if(v.i==v.Ncells-1 && v.right==NULL){
      const d& vVf=v.vVf;
      drhoim2+=-d_l()*vVf;
    }
    return drhoim2;
  }
  dmat Continuity::drhoip2(){
    dmat drhoip2=zero;
    if(v.i==0 && v.left==NULL){
      const d& vVf=v.vVf;
      drhoip2+=-d_r()*vVf;
    }
    // TRACE(10,"d_r():\n"<<d_r());
    return drhoip2;

  }
  Continuity::~Continuity(){}

} // Namespace tube

