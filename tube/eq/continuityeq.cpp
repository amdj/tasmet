#include "continuityeq.h"
#include "tubevertex.h"

#define CONT_VISCOSITY
#define CONT_SCALE (1.0)//(pow(v.gc.c0,2)) //(pow(v.gc.c0,2)/v.gc.omg)
#define CONT_SCALE0 (1.0)//(1.0/pow(v.gc.M,2))



namespace tube{
  vd Continuity::error(const TubeVertex& v) const {	// Current error in continuity equation
    // The default boundary implementation is an adiabatic no-slip wall.
    TRACE(6,"Continuity::Error()");
    vd error(v.gc->Ns,fillwith::zeros);
    error+=v.cWddt*v.gc->DDTfd*v.rho();
    error+=v.cWi*v.gc->fDFT*(v.rho.tdata()%v.U.tdata());
    if(v.i>0 || (v.i==0 && v.left!=NULL)){
      // Standard implementation of a no-slip (wall) boundary
      // condition
      vd rhoim1=v.left->rho.tdata();
      vd Uim1=v.left->U.tdata();
      error+=v.cWim1*v.gc->fDFT*(rhoim1%Uim1);
    }
    if(v.i<v.nCells-1 || (v.i==v.nCells-1 && v.right!=NULL) ){
      // Standard implementation of a no-slip (wall) boundary
      // condition
      vd rhoip1=v.right->rho.tdata();
      vd Uip1=v.right->U.tdata();
      error+=v.cWip1*v.gc->fDFT*(rhoip1%Uip1);
    }

    #ifdef CONT_VISCOSITY
    const d& vSf=v.lg.vSf;
    if(v.i>0 && v.i<v.nCells-1){
      error+=-d_r(v)*(v.right->rho() -v.rho())*vSf;
      error+= d_l(v)*(v.rho() -v.left->rho())  *vSf;
    }
    else if(v.i==0){		// First v
      error+=-d_r(v)*(v.right->right->rho()-v.right->rho())*vSf;
      error+=d_l(v)*(v.right->rho()-v.rho())*vSf;
    }
    else {			// Last v
      error+=-d_r(v)*(v.rho()-v.left->rho())*vSf;
      error+=d_l(v)*(v.left->rho()-v.left->left->rho())*vSf;
    }
    #endif
    // (Boundary) source term
    error+=v.csource();
    return error;
  }
  dmat Continuity::drhoi(const TubeVertex& v) const {
    TRACE(0,"Continuity::drhoi()");
    dmat drhoi=v.cWddt*v.gc->DDTfd;		// Initialize and add first term
    drhoi+=v.cWi*v.gc->fDFT*v.U.diagt()*v.gc->iDFT;

    // Artificial viscosity terms
    #ifdef CONT_VISCOSITY
    const d& vVf=v.lg.vVf;
    const d& vSf=v.lg.vSf;
    if(v.i>0 && v.i<v.nCells-1){
      drhoi+=(d_l(v)+d_r(v))*vSf;	// Middle vertex
    }
    else if(v.i==0)
      drhoi+=-d_l(v)*vSf;	// First vertex
    else		
      drhoi+=-d_r(v)*vSf;	// Last vertex
    #endif
    return drhoi;
  }
  dmat Continuity::dUi(const TubeVertex& v) const {
    TRACE(0,"Continuity::dUi()");
    dmat dUi=v.zero;
    dUi+=v.cWi*v.gc->fDFT*v.rho.diagt()*v.gc->iDFT;
    return dUi;
  }
  dmat Continuity::dUip1(const TubeVertex& v) const {
    TRACE(0,"Continuity::dUip1()");
    dmat dUip1=v.zero;
    if(v.right!=NULL)
      dUip1+=v.cWip1*v.gc->fDFT*v.right->rho.diagt()*v.gc->iDFT;
    return dUip1;
  }
  dmat Continuity::dUim1(const TubeVertex& v) const {
    TRACE(0,"Continuity::dUim1()");

    dmat dUim1=v.zero;
    if(v.left!=NULL)
      dUim1+=v.cWim1*v.gc->fDFT*v.left->rho.diagt()*v.gc->iDFT;
    return dUim1;
  }
  dmat Continuity::drhoip1(const TubeVertex& v) const {
    TRACE(0,"Continuity::drhoip1()");
    dmat drhoip1=v.zero;

    if(v.i<v.nCells-1 || v.right!=NULL)
      drhoip1=v.cWip1*v.gc->fDFT*v.right->U.diagt()*v.gc->iDFT;
    // Artificial viscosity terms
    #ifdef CONT_VISCOSITY    
    const d& vSf=v.lg.vSf;

    if(v.i>0 && v.i<v.nCells-1){
      drhoip1+=-d_r(v)*vSf;
    }
    else if(v.i==0){		// First vertex
      drhoip1+=(d_l(v)+d_r(v))*vSf;

    }
    #endif
    return drhoip1;
  }
  dmat Continuity::drhoim1(const TubeVertex& v) const {
    TRACE(0,"Continuity::drhoim1()");

    dmat drhoim1=v.zero;
    if(v.left!=NULL)
      drhoim1+=v.cWim1*v.gc->fDFT*v.left->U.diagt()*v.gc->iDFT;
    // Artificial viscosity terms
    #ifdef CONT_VISCOSITY
    const d& vSf=v.lg.vSf;
    if(v.i>0 && v.i<v.nCells-1){
      drhoim1+=-d_l(v)*vSf;
    }
    else if(v.i==v.nCells-1){		// Last vertex
      drhoim1+=(d_l(v)+d_r(v))*vSf;
    }
    #endif
    return drhoim1;
  }
  dmat Continuity::drhoim2(const TubeVertex& v) const {
    dmat drhoim2=v.zero;
    if(v.i==v.nCells-1 && v.right==NULL){
      const d& vSf=v.lg.vSf;
      drhoim2+=-d_l(v)*vSf;
    }
    return drhoim2;
  }
  dmat Continuity::drhoip2(const TubeVertex& v) const {
    dmat drhoip2=v.zero;
    if(v.i==0 && v.left==NULL){
      const d& vSf=v.lg.vSf;
      drhoip2+=-d_r(v)*vSf;
    }
    return drhoip2;
  }
} // Namespace tube

