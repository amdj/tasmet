#include "bccell.h"
#include "weightfactors.h"
#include "jacrow.h"
#include "state.h"

// #define STATE_SCALE (1)
#define STATE_SCALE (1/v.gc->p0())
#define iDFT (v.gc->iDFT())
#define fDFT (v.gc->fDFT())

// #define STATE_SCALE (1)
namespace duct{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;


  vd State::error()
   const {
    TRACE(15,"State::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    error+=v.p()();
    error(0)+=v.gc->p0();	       // Add p0 part

    const vd& rhot=v.rho().tdata();
    const vd& Tt=v.T().tdata();

    d Rs=v.gc->gas().Rs();
    error+=-fDFT*(Rs*rhot%Tt);
    
    return STATE_SCALE*error;
  }
  StateBc::StateBc(const BcCell& v):Equation(v),v(v){}
  vd StateBc::error()
   const {
    TRACE(15,"StateBc::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    error+=v.pbc()();
    error(0)+=v.gc->p0();	       // Add p0 part

    const vd& rhot=v.rhobc().tdata();
    const vd& Tt=v.Tbc().tdata();

    d Rs=v.gc->gas().Rs();
    error+=-fDFT*(Rs*rhot%Tt);
    
    return STATE_SCALE*error;
  }
  JacRow StateBc::jac() const{
    TRACE(15,"StateBc::jac()");
    JacRow jac(dofnr,5);
    TRACE(0,"StateBc, dofnr jac:"<< dofnr);
    jac+=JacCol(v.p(),STATE_SCALE*eye());

    const vd& rhot=v.rhobc().tdata();
    const vd& Tt=v.Tbc().tdata();

    d Rs=v.gc->gas().Rs();
    jac+=JacCol(v.rhobc(),-STATE_SCALE*Rs*fDFT*diagmat(Tt)*iDFT);
    jac+=JacCol(v.Tbc(),-STATE_SCALE*Rs*fDFT*diagmat(rhot)*iDFT);

    return jac;
  }
  JacRow State::jac() const{
    TRACE(15,"State::jac()");
    JacRow jac(dofnr,5);
    TRACE(0,"State, dofnr jac:"<< dofnr);
    jac+=JacCol(v.p(),STATE_SCALE*eye());

    const vd& rhot=v.rho().tdata();
    const vd& Tt=v.T().tdata();

    // vd rhomid=Wr*rhot+Wl*rhotL;
    // vd Tmid=Wr*Tt+Wl*TtL;
    d Rs=v.gc->gas().Rs();

    jac+=JacCol(v.rho(),-STATE_SCALE*Rs*fDFT*diagmat(Tt)*iDFT);
    jac+=JacCol(v.T(),-STATE_SCALE*Rs*fDFT*diagmat(rhot)*iDFT);

    return jac;
  }
  void State::show() const {
    cout << "----------------- State equation\n";
  }
  void StateBc::show() const {
    cout << "----------------- StateBc equation\n";
  }
} // namespace duct




