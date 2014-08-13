#include "momentumeq.h"
// #include "drag.h"
#include "tube.h"
#include "tubevertex.h"
// #define MOM_VISCOSITY

namespace tube{

  void Momentum::show() const{
    // cout << "--------- Showing momentum weight factors for i=" << v.i <<"\n" \
    // 	 << "v.mWddt      : "<<v.mWddt      <<"\n"					\
    // 	 << "v.mWuim1     : "<<v.mWuim1<<"\n"						\
    // 	 << "v.mWui       : "<<v.mWui      <<"\n"					\
    // 	 << "v.mWuip1     : "<<v.mWuip1<<"\n"						\
    // 	 << "v.mWpim1     : "<<v.mWpim1      <<"\n"					\
    // 	 << "v.mWpi       : "<<v.mWpi      <<"\n"					\
    // 	 << "v.mWpip1     : "<<v.mWpip1      <<"\n"					\
    //   ;      
  }
  void Momentum::init(const Tube& t)
  {
    TRACE(8,"Momentum::init(tube)");
    drag=&t.getDragResistance();
  }
  
  vd Momentum::error(const TubeVertex& v) const {		// Error in momentum equation

    TRACE(6,"Momentum::Error()");
    vd error(v.gc->Ns,fillwith::zeros);

    const vd& rhoti=v.rho.tdata();
    const vd& Uti=v.U.tdata();
    error+=v.mWddt*v.gc->DDTfd*v.gc->fDFT*(Uti%rhoti);
    error+=v.mWui*v.gc->fDFT*(rhoti%Uti%Uti);
    
    // Pressure terms    
    const vd& pi=v.p();
    error+=v.mWpi*pi;

    
    if(v.left!=NULL){
      const vd& rhotim1=v.left->rho.tdata();
      const vd& Utim1=v.left->U.tdata();
      error+=v.mWuim1*v.gc->fDFT*(rhotim1%Utim1%Utim1);
      // Pressure term
      const vd& pim1=v.left->p();
      error+=v.mWpim1*pim1;
    }
    if(v.right!=NULL){
      const vd& Utip1=v.right->U.tdata();
      const vd& rhotip1=v.right->rho.tdata();
      error+=v.mWuip1*v.gc->fDFT*(rhotip1%Utip1%Utip1);
      // Pressure term    
      const vd& pip1=v.right->p();
      error+=v.mWpip1*pip1;
    }

    // Artificial viscosity ter
    #ifdef MOM_VISCOSITY
    if(v.left!=NULL && v.right!=NULL){
      error+=-v.gc->rho0*d_r(v)*(v.right->U()-v.U());
      error+= v.gc->rho0*d_l(v)*(v.U()-v.left->U());
    }
    else if(v.left==NULL){
      error+=v.gc->rho0*d_r(v)*(v.right->U()-v.right->right->U());
      error+=v.gc->rho0*d_l(v)*(v.right->U()-v.U());
    }
    else if(v.right==NULL) {
      error+=v.gc->rho0*d_r(v)*(v.left->U()-v.U());
      error+=v.gc->rho0*d_l(v)*(v.left->U()-v.left->left->U());
    }
    #endif

    // Drag term
    assert(drag!=NULL);
    error+=v.mWddt*drag->drag(v);

    // (Boundary) source term
    error+=v.msource();
    return error;
  }
  dmat Momentum::dUi(const TubeVertex& v) const {
    TRACE(0,"Momentum::dUi()");
    dmat dUi=v.zero;
    // dUi+=vVf*tube.drag.dUi(i)/vSf;		       // Drag term
    dUi+=v.mWddt*v.gc->DDTfd*v.gc->fDFT*v.rho.diagt()*v.gc->iDFT; // Time-derivative term
    dUi+=2.0*v.mWui*v.gc->fDFT*(v.rho.diagt()*v.U.diagt())*v.gc->iDFT;
    // Artificial viscosity terms
    #ifdef MOM_VISCOSITY
    if(v.left!=NULL && v.right!=NULL){
      dUi+=v.gc->rho0*(d_l(v)+d_r(v));
    }
    else if(v.left==NULL)
      dUi+=-v.gc->rho0*d_l(v);
    else if(v.right==NULL)
      dUi+=-v.gc->rho0*d_r(v);
    #endif
    assert(drag!=NULL);
    dUi+=v.mWddt*drag->dUi(v);
    return dUi;
  }
  dmat Momentum::drhoi(const TubeVertex& v) const {
    TRACE(0,"Momentum::drhoi()");
    dmat drhoi=v.zero;
    drhoi+=v.mWddt*v.gc->DDTfd*v.gc->fDFT*v.U.diagt()*v.gc->iDFT;
    drhoi+=v.mWui*v.gc->fDFT*v.U.diagt()*v.U.diagt()*v.gc->iDFT;
    return drhoi;
  }
  dmat Momentum::dpi(const TubeVertex& v) const {
    TRACE(0,"Momentum::dpi()");
    dmat I(v.gc->Ns,v.gc->Ns,fillwith::eye);
    dmat dpi=v.zero;
    dpi+=v.mWpi*I;
    return dpi;
  }
  dmat Momentum::drhoim1(const TubeVertex& v) const {
    TRACE(0,"Momentum::drhoim1()");
    dmat drhoim1=v.zero;
    if(v.left!=NULL)
      drhoim1+=v.mWuim1*v.gc->fDFT*v.left->U.diagt()*v.left->U.diagt()*v.gc->iDFT;
    return drhoim1;
  }
  dmat Momentum::dUim1(const TubeVertex& v) const {
    TRACE(0,"Momentum::dUim1()");    // Todo: add this term!;
    dmat dUim1=v.zero;
    if(v.left!=NULL){
      dUim1+=2.0*v.mWuim1*v.gc->fDFT*v.left->rho.diagt()*v.left->U.diagt()*v.gc->iDFT;
    }
    #ifdef MOM_VISCOSITY
    // Artificial viscosity terms
    if(v.left!=NULL && v.right!=NULL){
      dUim1+=-v.gc->rho0*d_l(v);
    }
    else if(v.right==NULL)
      dUim1+=v.gc->rho0*(d_r(v)+d_l(v));
    #endif
    
    return dUim1;
  }
  dmat Momentum::dpim1(const TubeVertex& v) const {
    TRACE(0,"Momentum::dpim1()");
    dmat dpim1=v.zero;
    dmat I(v.gc->Ns,v.gc->Ns,fillwith::eye);
    if(v.left!=NULL)
      dpim1+=v.mWpim1*I;

    return dpim1;
  }
  dmat Momentum::drhoip1(const TubeVertex& v) const {
    TRACE(0,"Momentum::dhoip1()");    // Todo: add this term!;
    dmat drhoip1=v.zero;
    if(v.right!=NULL)
      drhoip1+=v.mWuip1*v.gc->fDFT*v.right->U.diagt()*v.right->U.diagt()*v.gc->iDFT;
    return drhoip1;
  }
  dmat Momentum::dUip1(const TubeVertex& v) const {
    TRACE(0,"Momentum::dUip1()"); // Todo: add this term!;
    dmat dUip1=v.zero;
    if(v.right!=NULL){
      dUip1+=2.0*v.mWuip1*v.gc->fDFT*v.right->rho.diagt()*v.right->U.diagt()*v.gc->iDFT;
    }
    #ifdef MOM_VISCOSITY
    // Artificial viscosity terms
    if(v.left!=NULL && v.right!=NULL){
      dUip1+=-v.gc->rho0*d_r(v);
    }
    else if(v.left==NULL)
      dUip1+=v.gc->rho0*(d_r(v)+d_l(v));
    #endif
    return dUip1;
  }
  dmat Momentum::dpip1(const TubeVertex& v) const {
    TRACE(0,"Momentum::dpip1()");
    dmat dpip1=v.zero;
    dmat I(v.gc->Ns,v.gc->Ns,fillwith::eye);
    if(v.right!=NULL)
      dpip1+=v.mWpip1*I;
    return dpip1;
  }
  dmat Momentum::dUip2(const TubeVertex& v) const {
    TRACE(0,"Momentum:dUip2()");
    // TRACE(50,"i:"<<i);
    dmat dUip2=v.zero;
    #ifdef MOM_VISCOSITY
    if(v.left==NULL)
      dUip2+=-v.gc->rho0*d_r(v);
    #endif
    return dUip2;
  }
  dmat Momentum::dUim2(const TubeVertex& v) const {
    dmat dUim2=v.zero;
    #ifdef MOM_VISCOSITY
    if(v.right==NULL )
      dUim2+=-v.gc->rho0*d_l(v);
    #endif
    return dUim2;
  }  

} // namespace tube
















