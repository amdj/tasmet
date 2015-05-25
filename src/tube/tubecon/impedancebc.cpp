// impedancebc.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "impedancebc.h"
#include "tube.h"
#include "bccell.h"
#include "jacobian.h"

#define Ns (gc->Ns())
#define eye (eye(Ns,Ns))

namespace tube {
  using namespace tasystem;
  
  ImpedanceBc::ImpedanceBc(us segnr,Pos position,const var& Z,d T0):
    TubeBc(segnr,position),
    Z(Z),
    T0(T0)
  {
    TRACE(15,"ImpedanceBc::ImpedanceBc()");
  }
  ImpedanceBc::ImpedanceBc(const ImpedanceBc& other,const TaSystem& sys):
    TubeBc(other,sys),
    Z(other.Z),
    T0(other.T0)
  {
    TRACE(15,"ImpedanceBc::ImpedanceBc(copy and init)");
    assert(gc);
    Z.setGc(*gc);
    setInit(true);
  }
  void ImpedanceBc::updateNf() {
    Z.updateNf();
  }
  void ImpedanceBc::setEqNrs(us firsteqnr) {
    TRACE(15,"void ImpedanceBc::setEqNrs()");
    TubeBc::setEqNrs(firsteqnr);
  } 
  vd ImpedanceBc::error() const {
    TRACE(15,"ImpedanceBc::error()");
    vd error(getNEqs(),fillwith::zeros);
    const BcCell& cell=t->bcCell(pos);


    // error.subvec(2*Ns,3*Ns-1)=-cell.mHbc()()
    // +fDFT*(cp(cell)*cell.mbc().tdata()%cell.Tbc().tdata());
    // VARTRACE(30,cell.Tbc()());
    // VARTRACE(30,cp(cell));
    error.subvec(2*Ns,3*Ns-1)=-cell.mHbc()();
    error.subvec(2*Ns,3*Ns-1)+=cell.extrapolateQuant(Varnr::mH);
    return error;
  }
  void ImpedanceBc::jac(Jacobian& jac) const{
    TRACE(8,"ImpedanceBc::jac()");

    const BcCell& cell=t->bcCell(pos);
    // Prescribed temperature Jacobian part

    // Momentum equation Jacobian
    // Prescribed enthalpy flow.
    JacRow enthalpy_extrapolated_jac(firsteqnr+2*Ns,3);

    enthalpy_extrapolated_jac+=cell.dExtrapolateQuant(Varnr::mH);
    enthalpy_extrapolated_jac+=JacCol(cell.mHbc(),-eye);
    jac+=enthalpy_extrapolated_jac;

  }
  void ImpedanceBc::show(us detailnr) const {
    TRACE(5,"ImpedanceBc::show()");
    checkInit();
    const char* side;
    if(pos==Pos::left)
      side="left";
    else
      side="right";
    cout << "ImpedanceBc boundary condition set at "<<side <<" side of segment "<< segnr<<".\n";
    if(detailnr>1){
      cout << "Prescribed impedance:" << Z.getcRes() << "\n";
      cout << "Prescribed fluid temperature:" << T0 << "\n";      
      // cout << "Prescribed solid temperature:" << prescribeTs.getVals()() << "\n";      
    }
  } // show

} // namespace tube

//////////////////////////////////////////////////////////////////////
