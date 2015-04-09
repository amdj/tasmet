// mHEq.cpp
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
  using std::ignore;

  namespace 
  {
    inline d cp(const Cell& c) {
      return c.gc->gas().cp(c.gc->T0());
    }
  } // namespace 

  void mHEq::init() {
    TRACE(15,"mHEq::init()");
    std::tie(WLl,WLr,ignore,ignore)=WeightFactors(v);
  }
  void mHEq::show() const{
    cout << "----------------- mHEq equation\n";
    cout << "WLl     :"<<WLl<<"\n";
    cout << "WLr     :"<<WLr<<"\n";
  }
  inline vd TLt_(const Cell& v,d WLl,d WLr){
    return WLl*v.TL().tdata()+WLr*v.T().tdata();
  }

  inline vd half_u_sq_td(const Cell& v){
    const vd& rhot=v.rho().tdata();
    const vd& mut=v.mu().tdata();
    d vSf=v.vSf;
    return 0.5*pow(mut/(rhot*vSf),2);
  }
  vd mHEq::error() const{
    TRACE(15,"mHEq::error()");
    assert(v.left());
    d cp0=cp(v);
    vd error=-v.mHL()();
    // Static enthalpy flow
    vd TLt=TLt_(v,WLl,WLr);
    error+=cp0*fDFT*(v.mL().tdata()%TLt);

    // Weighted average of kinetic energy flow
    // .5*u^2 left:

    // // u here
    if(v.i==1)
      WARN("Incomplete! Missing terms of kinetic energy. In Jacobian as well");
    error+=fDFT*(v.mL().tdata()%
                 (
                  WLl*half_u_sq_td(*v.left()) // From left side
                  +WLr*half_u_sq_td(v)        // From this side
                  )
                 );

    return error;
  }
  tasystem::JacRow mHEq::jac() const {
    TRACE(15,"mHEq::jac()");
    assert(v.left());
    JacRow jac(dofnr,4);
    d cp0=cp(v);

    // Static enthalpy terms
    const vd& mLt=v.mL().tdata();
    // TL time domain
    vd TLt=TLt_(v,WLl,WLr);

    jac+=JacCol(v.mHL(),-eye());
    jac+=JacCol(v.TL(),fDFT*diagmat(cp0*WLl*mLt)*iDFT);
    jac+=JacCol(v.T(),fDFT*diagmat(cp0*WLr*mLt)*iDFT);
    jac+=JacCol(v.mL(),fDFT*diagmat(cp0*TLt)*iDFT);

    // Kinetic energy terms
    vd halfu_sq_td=WLl*half_u_sq_td(*v.left()) // From left side
      +WLr*half_u_sq_td(v); // From this side

    jac+=JacCol(v.mL(),fDFT*diagmat(halfu_sq_td)*iDFT);

    return jac;
  }
  
  
} // namespace tube

//////////////////////////////////////////////////////////////////////
