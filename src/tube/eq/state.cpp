#include "tubevertex.h"
#include "weightfactors.h"
#include "jacobian.h"
#include "state.h"

#define STATE_SCALE (1/v.gc->p0)
#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)

// #define STATE_SCALE (1)
namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;


  // State class
  vd State::error()
   const {
    TRACE(6,"State::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    error+=v.pL()();
    error(0)+=v.gc->p0;	       // Add p0 part

    const vd& rhot=v.rho().tdata();
    const vd& Tt=v.T().tdata();
    d Rs=v.gc->gas.Rs();
    error+=fDFT*(WLi*Rs*rhot%Tt);
    error+=fDFT*(WLim1*Rs*v.rhoL().tdata()%v.TL().tdata());

    return STATE_SCALE*error;
  }
  void State::init() {
    const WeightFactors& w=v.weightFactors();
    if(v.left()){
      WLi=-w.wLr;
      WLim1=-w.wLl;
    }
    else{
      WLi=0;
      WLim1=-1;
    }
  }
  JacRow State::jac() const{
    TRACE(6,"State::jac()");
    JacRow jac(dofnr,5);
    TRACE(0,"State, dofnr jac:"<< dofnr);
    jac+=dpL();

    jac+=dTi();
    jac+=drhoi();

    jac+=dTim1();
    jac+=drhoim1();
    return jac;
  }
  JacCol State::dpL()
   const {
    TRACE(0,"State::dpL");
    return JacCol(v.pL(),STATE_SCALE*eye<dmat>(v.gc->Ns(),v.gc->Ns()));
  }
  JacCol State::dTi()
   const {
    TRACE(0,"State::dTi()");
    dmat rhotidiag=diagmat(v.rho().tdata());
    return JacCol(v.T(),WLi*STATE_SCALE*v.gc->gas.Rs()*fDFT*rhotidiag*iDFT);
  }
  JacCol State::drhoi()
   const {
    TRACE(0,"State::drhoi()");
    dmat Ttidiag=diagmat(v.T().tdata());
    return JacCol(v.rho(),WLi*STATE_SCALE*v.gc->gas.Rs()*fDFT*Ttidiag*iDFT);
  }
  JacCol State::dTim1()
   const {
    TRACE(0,"State::dTim1()");
    dmat rhotim1diag=diagmat(v.rhoL().tdata());
    return JacCol(v.TL(),WLim1*STATE_SCALE*v.gc->gas.Rs()*fDFT*rhotim1diag*iDFT);
  }
  JacCol State::drhoim1()
   const {
    TRACE(0,"State::drhoim1()");
    dmat Ttim1diag=diagmat(v.TL().tdata());
    return JacCol(v.rhoL(),WLim1*STATE_SCALE*v.gc->gas.Rs()*fDFT*Ttim1diag*iDFT);
  }

  // StateR class
  vd StateR::error()
   const {
    TRACE(6,"StateR::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    error+=v.pR()();
    error(0)+=v.gc->p0;	       // Add p0 part

    const vd& rhoRt=v.rhoR().tdata();
    const vd& TRt=v.TR().tdata();
    d Rs=v.gc->gas.Rs();
    error-=fDFT*(Rs*rhoRt%TRt);
    return STATE_SCALE*error;
  }
  void StateR::init(){
    TRACE(15,"StateR::init()");
  }
  JacRow StateR::jac() const{
    TRACE(6,"State::jac()");
    JacRow jac(dofnr,3);
    jac+=dpR();
    jac+=dTR();
    jac+=drhoR();
    return jac;
  }
  JacCol StateR::dpR()
   const {
    TRACE(0,"StateR::dpR");
    return JacCol(v.pR(),STATE_SCALE*eye<dmat>(v.gc->Ns(),v.gc->Ns()));
  }
  JacCol StateR::dTR()
   const {
    TRACE(0,"StateR::dTR()");
    dmat rhotidiag=diagmat(v.rho().tdata());
    return JacCol(v.T(),-STATE_SCALE*v.gc->gas.Rs()*fDFT*v.rhoR().diagt()*iDFT);
  }
  JacCol StateR::drhoR()  const {
    TRACE(0,"StateR::drhoip1()");
    return JacCol(v.rhoR(),-STATE_SCALE*v.gc->gas.Rs()*fDFT*v.TR().diagt()*iDFT);
  }

} // namespace tube




