// mu.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "mu.h"
#include "cell.h"
#include "jacrow.h"

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)
#define DDTfd (v.gc->DDTfd)

namespace tube {
  using tasystem::JacRow;
  using tasystem::JacCol;

  inline vd mu_td(const Cell& v) {
    return 0.5*(pow(v.mL().tdata(),2)+pow(v.mR().tdata(),2))/(v.vSf*v.rho().tdata());
  }
    
  vd MuEq::error() const{
    TRACE(15,"MuEq::error()");
    return -v.mu()()+fDFT*mu_td(v);
  }
  void MuEq::show() const{
    cout << "----------------- MuEq\n";
  }
  tasystem::JacRow MuEq::jac() const {
    TRACE(15,"MuEq::jac()");
    JacRow jac(dofnr,4);
    jac+=JacCol(v.mu(),-eye());

    const vd& rhot=v.rho().tdata();


    jac+=JacCol(v.rho(),-fDFT*diagmat(mu_td(v)/rhot)*iDFT);
    jac+=JacCol(v.mL(),fDFT*diagmat(0.5*v.mL().tdata()/(v.vSf*rhot))*iDFT);
    jac+=JacCol(v.mR(),fDFT*diagmat(0.5*v.mR().tdata()/(v.vSf*rhot))*iDFT);
    return jac;
  }
  
  
} // namespace tube

//////////////////////////////////////////////////////////////////////
