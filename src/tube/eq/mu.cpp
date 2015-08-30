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
    d Wfo=0;//v.gc->getWfo();
    return 0.5*((1-Wfo)*pow(v.ml().tdata(),2)+(1+Wfo)*pow(v.mr().tdata(),2))/(v.vSf*v.rho().tdata());
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

    d Wfo=0;//v.gc->getWfo();    
    jac+=JacCol(v.rho(),-fDFT*diagmat(mu_td(v)/rhot)*iDFT);
    jac+=JacCol(v.ml(),fDFT*diagmat(0.5*(1-Wfo)*v.ml().tdata()/(v.vSf*rhot))*iDFT);
    jac+=JacCol(v.mr(),fDFT*diagmat(0.5*(1+Wfo)*v.mr().tdata()/(v.vSf*rhot))*iDFT);
    return jac;
  }
  
  
} // namespace tube

//////////////////////////////////////////////////////////////////////
