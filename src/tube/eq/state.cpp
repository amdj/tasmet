#include "tubevertex.h"
#include "weightfactors.h"
#include "jacobian.h"
#include "state.h"

// #define STATE_SCALE (1)
#define STATE_SCALE (1/v.gc->p0)
#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)

// #define STATE_SCALE (1)
namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;


  // State class
  void State::init() {
    const WeightFactors& w=v.weightFactors();
    if(v.left()){
      Wl=w.wLl;
      Wr=w.wLr;
    }
    else{
      Wr=0;
      Wl=1;
    }
  }
  vd State::error()
   const {
    TRACE(6,"State::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    error+=v.pL()();
    error(0)+=v.gc->p0;	       // Add p0 part

    const vd& rhot=v.rho().tdata();
    const vd& Tt=v.T().tdata();
    const vd& rhotL=v.rhoL().tdata();
    const vd& TtL=v.TL().tdata();

    // vd rhomid=Wr*rhot+Wl*rhotL;
    // vd Tmid=Wr*Tt+Wl*TtL;
    d Rs=v.gc->gas.Rs();
    error+=-fDFT*(Wr*Rs*rhot%Tt);
    error+=-fDFT*(Wl*Rs*rhotL%TtL);
    
    return STATE_SCALE*error;
  }
  JacRow State::jac() const{
    TRACE(6,"State::jac()");
    JacRow jac(dofnr,5);
    TRACE(0,"State, dofnr jac:"<< dofnr);
    jac+=JacCol(v.pL(),STATE_SCALE*eye<dmat>(v.gc->Ns(),v.gc->Ns()));

    const vd& rhot=v.rho().tdata();
    const vd& Tt=v.T().tdata();
    const vd& rhotL=v.rhoL().tdata();
    const vd& TtL=v.TL().tdata();

    // vd rhomid=Wr*rhot+Wl*rhotL;
    // vd Tmid=Wr*Tt+Wl*TtL;
    d Rs=v.gc->gas.Rs();

    jac+=JacCol(v.rhoL(),-STATE_SCALE*Rs*Wl*fDFT*diagmat(TtL)*iDFT);
    jac+=JacCol(v.rho(),-STATE_SCALE*Rs*Wr*fDFT*diagmat(Tt)*iDFT);
    jac+=JacCol(v.TL(),-STATE_SCALE*Rs*Wl*fDFT*diagmat(rhotL)*iDFT);
    jac+=JacCol(v.T(),-STATE_SCALE*Rs*Wr*fDFT*diagmat(rhot)*iDFT);

    return jac;
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
    jac+= JacCol(v.pR(),STATE_SCALE*eye<dmat>(v.gc->Ns(),v.gc->Ns()));
    jac+=JacCol(v.T(),-STATE_SCALE*v.gc->gas.Rs()*fDFT*v.rhoR().diagt()*iDFT);
    jac+=JacCol(v.rhoR(),-STATE_SCALE*v.gc->gas.Rs()*fDFT*v.TR().diagt()*iDFT);
    return jac;
  }

} // namespace tube




