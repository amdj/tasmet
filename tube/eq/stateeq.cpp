#include "stateeq.h"
#include "tubevertex.h"
#include "jacobian.h"
#include "weightfactors.h"
#define STATE_SCALE (1/v.gc->p0)
// #define STATE_SCALE (1)
namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  // vd State::error()
  //  const {
  //   TRACE(6,"State::Error()");
  //   vd error(v.gc->Ns(),fillwith::zeros);
  //   // TRACE(-1,"State p0:"<<p0);
  //   error+=0.5*(v.pL()()+v.pR()());
  //   error(0)+=v.gc->p0;	       // Add p0 part
  //   // TRACE(-1,"state error:"<<error);    
  //   // TRACE(-1,"T0:"<<vertex.gc->gas.Rs()*fDFT()*(vertex.T.tdata()%vertex.rho().tdata()));    
  //   error+=-1.0*v.gc->gas.Rs()*v.gc->fDFT*(v.rho().tdata()%v.T().tdata());
  //   // TRACE(-1,"state error:"<<error);
  //   return STATE_SCALE*error;
  // }
  // JacRow State::jac() const{
  //   TRACE(6,"State::jac()");
  //   JacRow jac(dofnr,3);
  //   TRACE(0,"State, dofnr jac:"<< dofnr);    
  //   jac+=dpL();
  //   jac+=dpR();    
  //   jac+=dTi();
  //   jac+=drhoi();
  //   return jac;
  // }
  // JacCol State::dpL() const {
  //   TRACE(0,"State::dpi");
  //   return JacCol(v.pL(),0.5*STATE_SCALE*eye<dmat>(v.gc->Ns(),v.gc->Ns()));
  // }
  // JacCol State::dpR() const {
  //   TRACE(0,"State::dpi");
  //   return JacCol(v.pR(),0.5*STATE_SCALE*eye<dmat>(v.gc->Ns(),v.gc->Ns()));
  // }
  // JacCol State::dTi()
  //  const {
  //   TRACE(10,"State::dTi()");
  //   dmat rhotidiag=diagmat(v.rho().tdata());
  //   // TRACE(15,"dTi value:"<< -1.0*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*rhotidiag*v.gc->iDFT);
  //   return JacCol(v.T(),-1.0*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*rhotidiag*v.gc->iDFT);
  // }
  // JacCol State::drhoi()
  //  const {
  //   TRACE(10,"State::drhoi()");
  //   dmat Ttidiag=diagmat(v.T().tdata());
  //   // TRACE(15,"dTi value:"<< -1.0*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*Ttidiag*v.gc->iDFT);
  //   return JacCol(v.rho(),-1.0*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*Ttidiag*v.gc->iDFT);
  // }

  vd StateL::error()
   const {
    TRACE(6,"StateL::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    // TRACE(-1,"StateL p0:"<<p0);
    error+=v.pL()();
    error(0)+=v.gc->p0;	       // Add p0 part
    vd rhoTti=v.rho().tdata()%v.T().tdata();
    error+=v.gc->gas.Rs()*(v.gc->fDFT*(WLi*rhoTti));
    if(v.left())  
      error+=v.gc->gas.Rs()*(v.gc->fDFT*(WLim1*v.rhoL().tdata()%v.TL().tdata()));
    else
      error+=v.gc->gas.Rs()*(v.gc->fDFT*(WLip1*v.right()->rho().tdata()%v.TR().tdata()));      
    return STATE_SCALE*error;
  }
  void StateL::init() {

    if(v.left()){
      WLi=-w.wLr;
      WLim1=-w.wLl;
      WLip1=0;
    }
    else{
      WLi=-w.wL0;
      WLim1=0;
      WLip1=-w.wL1;
    }

  }

  JacRow StateL::jac() const{
    TRACE(6,"StateL::jac()");
    JacRow jac(dofnr,5);
    TRACE(0,"StateL, dofnr jac:"<< dofnr);
    jac+=dpL();
    jac+=dTi();
    jac+=drhoi();
    if(v.left()){
      jac+=dTim1();
      jac+=drhoim1();
    }
    else{
      jac+=dTip1();
      jac+=drhoip1();
    }
    return jac;
  }
  JacCol StateL::dpL()
   const {
    TRACE(0,"StateL::dpi");
    return JacCol(v.pL(),STATE_SCALE*eye<dmat>(v.gc->Ns(),v.gc->Ns()));
  }

  JacCol StateL::dTi()
   const {
    TRACE(0,"StateL::dTi()");
    dmat rhotidiag=diagmat(v.rho().tdata());
    return  JacCol(v.T(),WLi*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*rhotidiag*v.gc->iDFT);
  }
  JacCol StateL::drhoi()
   const {
    TRACE(0,"StateL::drhoi()");
    dmat Ttidiag=diagmat(v.T().tdata());
    return JacCol(v.rho(),WLi*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*Ttidiag*v.gc->iDFT);
  }
  JacCol StateL::dTim1()
   const {
    TRACE(0,"StateL::dTim1()");
    dmat rhotim1diag=diagmat(v.rhoL().tdata());
    return JacCol(v.TL(),WLim1*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*rhotim1diag*v.gc->iDFT);
  }
  JacCol StateL::drhoim1()
   const {
    TRACE(0,"StateL::drhoim1()");
    dmat Ttim1diag=diagmat(v.TL().tdata());
    return JacCol(v.rhoL(),WLim1*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*Ttim1diag*v.gc->iDFT);
  }
  JacCol StateL::dTip1()
   const {
    TRACE(0,"StateL::dTip1()");
    dmat rhotip1diag=diagmat(v.right()->rho().tdata());
    return JacCol(v.TR(),WLip1*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*rhotip1diag*v.gc->iDFT);
  }
  JacCol StateL::drhoip1()
   const {
    TRACE(0,"StateL::drhoip1()");
    dmat Ttip1diag=diagmat(v.TR().tdata());
    return JacCol(v.rhoR(),WLip1*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*Ttip1diag*v.gc->iDFT);
  }
    
  // *********************************************** Now for pR,
  // *********************************************** generally a last
  // *********************************************** node

  //   vd StateR::error()
  //  const {
  //   TRACE(6,"StateR::Error()");
  //   vd error(v.gc->Ns(),fillwith::zeros);
  //   // TRACE(-1,"StateL p0:"<<p0);
  //   error+=v.pR()();
  //   error(0)+=v.gc->p0;	       // Add p0 part
  //   vd rhoTti=v.rho().tdata()%v.T().tdata();
  //   error+=-v.gc->gas.Rs()*(v.gc->fDFT*(wRNm1*rhoTti));
  //   if(v.left())  
  //     error+=-v.gc->gas.Rs()*(v.gc->fDFT*(wRNm2*v.rhoL().tdata()%v.T.tdataL()));
  //   else{
  //     WARN("StateR equation not implemented correctly");
  //     exit(1);
  //   }
  //   return STATE_SCALE*error;
  // }
  // JacRow StateR::jac() const{
  //   TRACE(6,"StateR::jac()");
  //   JacRow jac(dofnr,5);
  //   TRACE(0,"StateR, dofnr jac:"<< dofnr);
  //   jac+=dpR();
  //   jac+=dTi();
  //   jac+=drhoi();
  //   if(v.left){
  //     jac+=dTim1();
  //     jac+=drhoim1();
  //   }
  //   else{
  //     WARN("StateR equation not implemented correctly");
  //     exit(1);
  //   }
  //   return jac;
  // }
  // JacCol StateR::dpR()
  //  const {
  //   TRACE(10,"StateR::dpR");
  //   return JacCol(v.pR(),STATE_SCALE*eye<dmat>(v.gc->Ns(),v.gc->Ns()));
  // }
  // JacCol StateR::dTi()
  //  const {
  //   TRACE(10,"StateR::dTi()");
  //   dmat rhotidiag=diagmat(v.rho().tdata());
  //   return JacCol(v.T,-v.wRNm1*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*rhotidiag*v.gc->iDFT);
  // }    
  // JacCol StateR::drhoi()
  //  const {
  //   TRACE(10,"StateR::drhoi()");
  //   dmat Ttidiag=diagmat(v.T()tdata());
  //   return JacCol(v.rho,-v.wRNm1*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*Ttidiag*v.gc->iDFT);
  // }

  // JacCol StateR::dTim1()
  //  const {
  //   TRACE(10,"StateR::dTim1()");
  //   dmat rhotim1diag=diagmat(v.rhoL().tdata());
  //   return JacCol(v.TL(),-v.wRNm2*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*rhotim1diag*v.gc->iDFT);
  // }
  // JacCol StateR::drhoim1()
  //  const {
  //   TRACE(10,"StateR::drhoim1()");
  //   dmat Ttim1diag=diagmat(v.left()->T.tdata());
  //   return JacCol(v.rhoL(),-v.wRNm2*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*Ttim1diag*v.gc->iDFT);
  // }    
} // namespace tube




