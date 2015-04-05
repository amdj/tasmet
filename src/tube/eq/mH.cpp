// mH.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "mH.h"
#include "weightfactors.h"
#include "cell.h"
#include "jacrow.h"

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)
#define DDTfd (v.gc->DDTfd)
#define Ns (v.gc->Ns())

namespace tube {
  using tasystem::JacRow;
  using tasystem::JacCol;

  namespace 
  {
    inline d cp(const Cell& c) {
      return c.gc->gas().cp(c.gc->T0());
    }
  } // namespace 

  void mH::init() {
    TRACE(15,"mH::init()");
    std::tie(WLl,WLr,WRl,WRr)=WeightFactors(v);
  }
  void mH::show() const{
    cout << "------------- Full energy equation\n";
    cout << "WLl     :"<<WLl<<"\n";
    cout << "WLr     :"<<WLr<<"\n";
    cout << "WRl     :"<<WRl<<"\n";
    cout << "WRr     :"<<WRr<<"\n";
  }
  vd mH::error() const{
    TRACE(15,"mH::error()");
    assert(v.left());
    d cp0=cp(v);
    vd error=-v.mHL()();
    // Static enthalpy flow
    error+=cp0*fDFT*(v.mL().tdata()%(WLl*v.TL().tdata()+WLr*v.T().tdata()));

    // Weighted average of kinetic energy flow
    // 0.5*u^2 left:
    const vd& rhotL=v.left()->rho().tdata();
    const vd& muLt=v.left()->mu().tdata();
    d vSfL=v.left()->vSf;
    vd half_u_sq_l_td=0.5*muLt/(rhotL*vSfL);

    // 0.5*u^2 here
    const vd& rhot=v.rho().tdata();
    const vd& mut=v.mu().tdata();
    d vSf=v.vSf;
    vd half_u_sq_td=0.5*mut/(rhot*vSf);
    WARN("Incomplete! Missing terms of kinetic energy. In Jacobian as well");
    // error+=fDFT*(v.mL().tdata()%
                 // (WLl*half_u_sq_l_td+WLr*half_u_sq_td));

    return error;
  }
  tasystem::JacRow mH::jac() const {
    TRACE(15,"mH::jac()");
    JacRow jac(dofnr,4);
    d cp0=cp(v);
    const vd& mLt=v.mL().tdata();
    jac+=JacCol(v.mHL(),-eye());

    // TL time domain
    vd TLt=WLl*v.TL().tdata()+WLr*v.T().tdata();

    jac+=JacCol(v.TL(),fDFT*diagmat(cp0*WLl*mLt)*iDFT);
    jac+=JacCol(v.T(),fDFT*diagmat(cp0*WLr*mLt)*iDFT);
    jac+=JacCol(v.mL(),fDFT*diagmat(cp0*TLt)*iDFT);
    
    return jac;
  }
  
  
} // namespace tube

//////////////////////////////////////////////////////////////////////
