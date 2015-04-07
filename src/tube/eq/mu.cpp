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

  vd MuEq::error() const{
    TRACE(15,"MuEq::error()");
    return -v.mu()()+
      fDFT*(pow(0.5*(v.mL().tdata()+v.mR().tdata()),2)/(v.rho().tdata()*v.vSf));
  }
  void MuEq::show() const{
    cout << "----------------- MuEq\n";
  }
  tasystem::JacRow MuEq::jac() const {
    TRACE(15,"MuEq::jac()");
    JacRow jac(dofnr,4);
    jac+=JacCol(v.mu(),-eye());

    vd two_m_i_td=(v.mL().tdata()+v.mR().tdata());
    vd m_i_td_sq=pow(0.5*two_m_i_td,2);
    const vd& rhot=v.rho().tdata();
    jac+=JacCol(v.rho(),-fDFT*diagmat(m_i_td_sq/(v.vSf*pow(rhot,2)))*iDFT);
    jac+=JacCol(v.mL(),fDFT*diagmat(two_m_i_td/(v.vSf*rhot))*iDFT);
    jac+=JacCol(v.mR(),fDFT*diagmat(two_m_i_td/(v.vSf*rhot))*iDFT);
    return jac;
  }
  
  
} // namespace tube

//////////////////////////////////////////////////////////////////////
