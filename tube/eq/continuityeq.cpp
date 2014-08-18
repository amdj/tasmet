#include "continuityeq.h"
#include "tubevertex.h"

#include "artvisco.h"


namespace tube{
  void Continuity::show() const {
    cout << "Continuity equation\n";
    #ifdef CONT_VISCOSITY
    cout << "Artificial viscosity turned ON for continuity equation\n";
    #else
    cout << "Artificial viscosity turned OFF for continuity equation\n";
    #endif
    
  }
  void Continuity::init(const Tube& t){
    TRACE(8,"Continuity::init(tube)");
  }
    
  vd Continuity::error(const TubeVertex& v) const {	// Current error in continuity equation
    // The default boundary implementation is an adiabatic no-slip wall.
    TRACE(6,"Continuity::Error()");
    vd error(v.gc->Ns,fillwith::zeros);
    TRACE(10,"Ns:"<<v.gc->Ns);
    error+=v.cWddt*v.gc->DDTfd*v.rho();
    error+=v.cWi*v.gc->fDFT*(v.rho.tdata()%v.U.tdata());
    if(v.left!=NULL){
      // Standard implementation of a no-slip (wall) boundary
      // condition
      const vd& rhoim1=v.left->rho.tdata();
      const vd& Uim1=v.left->U.tdata();
      error+=v.cWim1*v.gc->fDFT*(rhoim1%Uim1);
    }
    // TRACE(10,"Right:,"<<v.right);    
    if(v.right!=NULL ){
      // Standard implementation of a no-slip (wall) boundary
      // condition
      const vd& rhoip1=v.right->rho.tdata();
      const vd& Uip1=v.right->U.tdata();
      error+=v.cWip1*v.gc->fDFT*(rhoip1%Uip1);
    }
    #ifdef CONT_VISCOSITY
    if(v.left!=NULL && v.right!=NULL){
      error+=d_l(v)*(v.cWart2*v.rho()+v.cWart1*v.left->rho() );
      error+=d_r(v)*(v.cWart3*v.rho()+v.cWart4*v.right->rho() );
    }
    #endif
    // (Boundary) source term
    error+=v.csource();
    TRACE(10,"SFSG");    
    return error;
  }
  dmat Continuity::drhoi(const TubeVertex& v) const {
    TRACE(0,"Continuity::drhoi()");
    dmat drhoi=v.cWddt*v.gc->DDTfd;		// Initialize and add first term
    drhoi+=v.cWi*v.gc->fDFT*v.U.diagt()*v.gc->iDFT;

    // Artificial viscosity terms
    #ifdef CONT_VISCOSITY
    if(v.left!=NULL && v.right!=NULL){
      drhoi+=(d_l(v)*v.cWart2+d_r(v)*v.cWart3);	// Middle vertex
    }
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
    if(v.right!=NULL)
      drhoip1=v.cWip1*v.gc->fDFT*v.right->U.diagt()*v.gc->iDFT;
    // Artificial viscosity terms
    #ifdef CONT_VISCOSITY
    if(v.left!=NULL && v.right!=NULL)
      drhoip1+=d_r(v)*v.cWart4;
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
    if(v.left!=NULL && v.right!=NULL)
      drhoim1+=d_l(v)*v.cWart1;
    #endif
    return drhoim1;
  }
  dmat Continuity::drhoim2(const TubeVertex& v) const {
    dmat drhoim2=v.zero;
    #ifdef CONT_VISCOSITY

    #endif    
    return drhoim2;
  }
  dmat Continuity::drhoip2(const TubeVertex& v) const {
    dmat drhoip2=v.zero;
    #ifdef CONT_VISCOSITY

    #endif
    return drhoip2;
  }
} // Namespace tube

