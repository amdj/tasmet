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
    const d& vSf=v.vSf;
    if(v.i>0 && v.i<v.Ncells-1){
      error+=-D_r()*(v.right->rho() -v.rho())*vSf;
      error+= D_l()*(v.rho() -v.left->rho())  *vSf;
    }
    else if(v.i==0){		// First v
      error+=-D_r()*(v.right->right->rho()-v.right->rho())*vSf;
      error+=D_l()*(v.right->rho()-v.rho())*vSf;
    }
    else {			// Last v
      error+=-D_r()*(v.rho()-v.left->rho())*vSf;
      error+=D_l()*(v.left->rho()-v.left->left->rho())*vSf;
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
      drhoi+=(D_l()+D_r())*vSf;	// Middle vertex
    }
    else if(v.i==0)
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
    const d& vSf=v.vSf;
    if(v.i>0 && v.i<v.Ncells-1){
      drhoip1+=-D_r()*vSf;
    }
    else if(v.i==0){		// First vertex
      drhoip1+=(D_l()+D_r())*vSf;
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
    const d& vSf=v.vSf;
    if(v.i>0 && v.i<v.Ncells-1){
      drhoim1+=-D_l()*vSf;
    }
    else if(v.i==v.Ncells-1){		// Last vertex
      drhoim1+=(D_l()+D_r())*vSf;
    }
    #endif

    drhoim1.row(0)*=CONT_SCALE0;
    return CONT_SCALE*drhoim1;
  }
  dmat Continuity::drhoim2(){
    dmat drhoim2=zero;
    if(v.i==v.Ncells-1 && v.right==NULL){
      const d& vSf=v.vSf;
      drhoim2+=-D_l()*vSf;
    }
    return drhoim2;
  }
  dmat Continuity::drhoip2(){
    dmat drhoip2=zero;
    if(v.i==0 && v.left==NULL){
      const d& vSf=v.vSf;
      drhoip2+=-D_r()*vSf;
    }
    // TRACE(10,"D_r():\n"<<D_r());
    return drhoip2;

  }
  Continuity::~Continuity(){}

} // Namespace tube

