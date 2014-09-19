// #define TRACERPLUS 20
#include "momentumeq.h"
#include "tube.h"
#include "tubevertex.h"

#include "artvisco.h"

namespace tube{

  void Momentum::show() const{
    cout << "Momentum equation\n";
    #ifdef MOM_VISCOSITY
    cout << "Artificial viscosity turned ON for momentum equation\n";
    #else
    cout << "Artificial viscosity turned OFF for momentum equation\n";
    #endif
    
  }
  void Momentum::init(const Tube& t)
  {
    TRACE(8,"Momentum::init(tube)");
    drag=&t.getDragResistance();
  }
  
  vd Momentum::error(const TubeVertex& v) const {		// Error in momentum equation

    TRACE(6,"Momentum::Error()");
    vd error(v.gc->Ns,fillwith::zeros);
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;
    const vd& rhoti=v.rho.tdata();
    const vd& Uti=v.U.tdata();
    error+=v.mWddt*v.gc->DDTfd*v.gc->fDFT*(Uti%rhoti);
    error+=v.mWui*v.gc->fDFT*(rhoti%Uti%Uti);
    
    // Pressure terms    
    error+=v.mWpL*v.pL()();
    error+=v.mWpR*v.pR()();    

    if(v.left!=NULL){
      const vd& rhotim1=v.left->rho.tdata();
      const vd& Utim1=v.left->U.tdata();
      error+=v.mWuim1*v.gc->fDFT*(rhotim1%Utim1%Utim1);
      // Pressure term
    }
    if(v.right!=NULL){
      const vd& Utip1=v.right->U.tdata();
      const vd& rhotip1=v.right->rho.tdata();
      error+=v.mWuip1*v.gc->fDFT*(rhotip1%Utip1%Utip1);
    }

    // Artificial viscosity ter
    #ifdef MOM_VISCOSITY
    if(v.left!=NULL && v.right!=NULL){
      error+=d_l(v)*fDFT*(v.mWart2*v.rho.tdata()%v.U.tdata()+v.mWart1*v.left->rho.tdata()%v.left->U.tdata());
      error+=d_r(v)*fDFT*(v.mWart3*v.rho.tdata()%v.U.tdata()+v.mWart4*v.right->rho.tdata()%v.right->U.tdata());
    }
    #endif

    // Drag term
    assert(drag!=NULL);
    error+=v.mWddt*drag->drag(v);

    // (Boundary) source term
    error+=v.msource();
    return error;
  }
  JacRow Momentum::jac(const TubeVertex& v) const {
    TRACE(6,"Momentum::jac()");
    JacRow jac(dofnr,9);
    TRACE(0,"Momentum, dofnr jac:"<< dofnr);
    jac+=drhoi(v);
    jac+=dUi(v);
    jac+=dpL(v);
    jac+=dpR(v);
    if(v.left){
      jac+=drhoim1(v);
      jac+=dUim1(v);
    }
    if(v.right){
      jac+=drhoip1(v);
      jac+=dUip1(v);
    }
    return jac;
  }
  void Momentum::domg(const TubeVertex& v,vd & domg_) const {
    TRACE(0,"Momentum::domg()");
    // Possibly later adding drag->domg();
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const us& Ns=v.gc->Ns;
    domg_.subvec(dofnr,dofnr+Ns-1)=v.gc->DDTfd*v.mWddt*fDFT*(v.rho.tdata()%v.U.tdata())/v.gc->getomg();
    TRACE(0,"Momentum::domg() done");
  }
  JacCol Momentum::dUi(const TubeVertex& v) const {
    TRACE(0,"Momentum::dUi()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;

    JacCol dUi(v.U);
    // dUi+=vVf*tube.drag.dUi(i)/vSf;		       // Drag term
    dUi+=v.mWddt*v.gc->DDTfd*v.gc->fDFT*v.rho.diagt()*v.gc->iDFT; // Time-derivative term
    dUi+=2.0*v.mWui*v.gc->fDFT*(v.rho.diagt()*v.U.diagt())*v.gc->iDFT;
    // Artificial viscosity terms
    #ifdef MOM_VISCOSITY
    if(v.left!=NULL && v.right!=NULL)
      dUi+=(v.mWart2*d_l(v)+v.mWart2*d_r(v))*fDFT*v.rho.diagt()*iDFT;
    #endif
    assert(drag!=NULL);
    dUi+=v.mWddt*drag->dUi(v);
    return dUi;
  }
  JacCol Momentum::drhoi(const TubeVertex& v) const {
    TRACE(0,"Momentum::drhoi()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;
    JacCol drhoi(v.rho);
    drhoi+=v.mWddt*v.gc->DDTfd*v.gc->fDFT*v.U.diagt()*v.gc->iDFT;
    drhoi+=v.mWui*v.gc->fDFT*v.U.diagt()*v.U.diagt()*v.gc->iDFT;
    #ifdef MOM_VISCOSITY
    if(v.left!=NULL && v.right!=NULL)
      drhoi+=(v.mWart2*d_l(v)+v.mWart3*d_r(v))*fDFT*v.U.diagt()*iDFT;
    #endif

    return drhoi;
  }
  JacCol Momentum::drhoim1(const TubeVertex& v) const {
    TRACE(0,"Momentum::drhoim1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;
    JacCol drhoim1(v.left->rho);
    drhoim1+=v.mWuim1*v.gc->fDFT*v.left->U.diagt()*v.left->U.diagt()*v.gc->iDFT;
    #ifdef MOM_VISCOSITY
    if(v.left!=NULL && v.right!=NULL)
      drhoim1+=v.mWart1*d_l(v)*fDFT*v.left->U.diagt()*iDFT;
    #endif

    return drhoim1;
  }
  JacCol Momentum::dUim1(const TubeVertex& v) const {
    TRACE(0,"Momentum::dUim1()");    // Todo: add this term!;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;

    JacCol dUim1(v.left->U);
    dUim1+=2.0*v.mWuim1*v.gc->fDFT*v.left->rho.diagt()*v.left->U.diagt()*v.gc->iDFT;
    // Artificial viscosity terms
    #ifdef MOM_VISCOSITY
    if(v.left!=NULL && v.right!=NULL){
      dUim1+=v.mWart1*d_l(v)*fDFT*v.left->rho.diagt()*iDFT;
    }
    #endif
    
    return dUim1;
  }
  JacCol Momentum::drhoip1(const TubeVertex& v) const {
    TRACE(0,"Momentum::dhoip1()");    // Todo: add this term!;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;
    JacCol drhoip1(v.right->rho);
    drhoip1+=v.mWuip1*v.gc->fDFT*v.right->U.diagt()*v.right->U.diagt()*v.gc->iDFT;
    #ifdef MOM_VISCOSITY
    // Artificial viscosity terms
    if(v.left!=NULL && v.right!=NULL){
      drhoip1+=v.mWart4*d_r(v)*fDFT*v.right->U.diagt()*iDFT;
    }
    #endif

    return drhoip1;
  }
  JacCol Momentum::dUip1(const TubeVertex& v) const {
    TRACE(0,"Momentum::dUip1()"); // Todo: add this term!;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;

    JacCol dUip1(v.right->U);
    dUip1+=2.0*v.mWuip1*v.gc->fDFT*v.right->rho.diagt()*v.right->U.diagt()*v.gc->iDFT;
    #ifdef MOM_VISCOSITY
    // Artificial viscosity terms
    if(v.left!=NULL && v.right!=NULL){
      dUip1+=v.mWart4*d_r(v)*fDFT*v.right->rho.diagt()*iDFT;
    }
    #endif
    return dUip1;
  }
  JacCol Momentum::dpR(const TubeVertex& v) const {
    TRACE(0,"Momentum::dpR()");

    dmat I(v.gc->Ns,v.gc->Ns,fillwith::eye);
    JacCol dpR(v.pR(),v.mWpR*I);
    return dpR;
  }
  JacCol Momentum::dpL(const TubeVertex& v) const {
    TRACE(0,"Momentum::dpR()");

    dmat I(v.gc->Ns,v.gc->Ns,fillwith::eye);
    JacCol dpL(v.pL(),v.mWpL*I);
    return dpL;
  }
  // JacCol Momentum::dUip2(const TubeVertex& v) const {
  //   TRACE(0,"Momentum:dUip2()");
  //   // TRACE(50,"i:"<<i);
  //   dmat dUip2=v.zero;
  //   #ifdef MOM_VISCOSITY
  //   // if(v.left==NULL)
  //     // dUip2+=-v.gc->rho0*d_l(v);
  //   #endif
  //   return dUip2;
  // }
  // JacCol Momentum::dUim2(const TubeVertex& v) const {
  //   dmat dUim2=v.zero;
  //   #ifdef MOM_VISCOSITY
  //   // if(v.right==NULL )
  //     // dUim2+=-v.gc->rho0*d_r(v);
  //   #endif
  //   return dUim2;
  // }  

} // namespace tube
