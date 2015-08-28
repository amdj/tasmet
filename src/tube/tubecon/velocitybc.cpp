// velocitybc.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "velocitybc.h"

#include "weightfactors.h"

#include "tube.h"
#include "bccell.h"
#include "jacobian.h"
#include "pycallback.h"
#include "Python.h"

#define fDFT (gc->fDFT)
#define iDFT (gc->iDFT)
#define DDTfd (gc->DDTfd)

#define Ns (gc->Ns())
#define Nf (gc->Nf())
#define eye (eye(Ns,Ns))

namespace tube {
  using namespace tasystem;

  VelocityBc::VelocityBc(const string& segid,Pos position,const var& up,d T0,bool arbitrateMass):
    TubeBc(segid,position),
    arbitrateMass(arbitrateMass),
    T0(T0),
    u_p(up)    
  {
    TRACE(15,"VelocityBc::VelocityBc()");

  }
  VelocityBc::VelocityBc(const VelocityBc& other,const TaSystem& sys):
    TubeBc(other,sys),
    arbitrateMass(other.arbitrateMass),
    T0(other.T0),
    u_p(other.u_p)
  {
    TRACE(15,"VelocityBc::VelocityBc(copy and init)");
    assert(gc);
    u_p.setGc(*gc);

    setInit(true);
  }
  VelocityBc::~VelocityBc(){}
  int VelocityBc::arbitrateMassEq() const {
    TRACE(15,"AdiabaticWall::arbitrateMassEq()");
    if(arbitrateMass)
      return firsteqnr;
    else
      return -1;
  }

  void VelocityBc::updateNf() {
    u_p.updateNf();
  }
  void VelocityBc::setEqNrs(us firsteqnr) {
    TRACE(15,"void VelocityBc::setEqNrs()");
    VARTRACE(10,firsteqnr);
    TubeBc::setEqNrs(firsteqnr);
  } 

  vd VelocityBc::error() const {
    TRACE(15,"VelocityBc::error()");
    vd error(getNEqs());
    const BcCell& cell=t->bcCell(pos);

    // Velocity error
    d Sf=pos==Pos::left?cell.SfL:cell.SfR;
    error.subvec(0,Ns-1)=fDFT*(cell.mbc().tdata()/(Sf*iDFT*cell.extrapolateQuant(Varnr::rho)))
      -u_p();

    // Isentropic state eq. error
    const vd p=cell.extrapolateQuant(Varnr::p);
    const var& Tbc=cell.Tbc();
    const d p0=gc->p0();
    const d gamma=gc->gas().gamma(T0);
    // Adiabatic compression/expansion
    error.subvec(Ns,2*Ns-1)=p-fDFT*(p0*(pow(Tbc.tdata()/T0,gamma/(gamma-1))-1));

    // Extrapolated enthalpy flow
    error.subvec(2*Ns,3*Ns-1)=-cell.mHbc()();
    error.subvec(2*Ns,3*Ns-1)+=cell.extrapolateQuant(Varnr::mH);
    return error;
  }
  void VelocityBc::jac(Jacobian& jac) const{
    TRACE(8,"VelocityBc::jac()");
    VARTRACE(10,firsteqnr);
    const BcCell& cell=t->bcCell(pos);

    // Velocity Jacobian
    // VARTRACE(25,cell.extrapolateQuant(Varnr::rho));

    JacRow impjac(firsteqnr,5); // 5=(2*p,2*rho,1*mbc)
    d Sf=pos==Pos::left?cell.SfL:cell.SfR;
    const vd rhobct=iDFT*cell.extrapolateQuant(Varnr::rho);
    impjac+=JacCol(cell.mbc(),fDFT*diagmat(1/(Sf*rhobct))*iDFT);
    d w1,w2; std::tie(w1,w2)=BcWeightFactorsV(cell);
    if(pos==Pos::left) {
      impjac+=JacCol(cell.rho(),-w1*fDFT*diagmat(cell.mbc().tdata()/(Sf*pow(rhobct,2)))*iDFT);      
      impjac+=JacCol(cell.rhoR(),-w2*fDFT*diagmat(cell.mbc().tdata()/(Sf*pow(rhobct,2)))*iDFT);      
    } else {
      impjac+=JacCol(cell.rho(),-w1*fDFT*diagmat(cell.mbc().tdata()/(Sf*pow(rhobct,2)))*iDFT);      
      impjac+=JacCol(cell.rhoL(),-w2*fDFT*diagmat(cell.mbc().tdata()/(Sf*pow(rhobct,2)))*iDFT);      
    }

    jac+=impjac;


    // Adiabatic pressure-temperature Jacobian:
    JacRow IsentropicJac(firsteqnr+Ns,2);
    const var& Tbc=cell.Tbc();

    // Pressure part
    IsentropicJac+=cell.dExtrapolateQuant(Varnr::p);
    const d p0=gc->p0();
    const d gamma=gc->gas().gamma(T0);

    // Temperature part
    d fac1=p0*gamma/(T0*(gamma-1));
    IsentropicJac+=JacCol(cell.Tbc(),-fDFT*diagmat(fac1*pow(Tbc.tdata()/T0,1/(gamma-1)))*iDFT);
    jac+=IsentropicJac;

    // Prescribed enthalpy flow.
    JacRow enthalpy_extrapolated_jac(firsteqnr+2*Ns,3);
    enthalpy_extrapolated_jac+=cell.dExtrapolateQuant(Varnr::mH);
    enthalpy_extrapolated_jac+=JacCol(cell.mHbc(),-eye);
    jac+=enthalpy_extrapolated_jac;

  }
  void VelocityBc::show(us detailnr) const {
    TRACE(5,"VelocityBc::show()");
    checkInit();
    const char* side=posWord(pos);

    cout << "VelocityBc boundary condition set at "<<side <<" side of segment "<< segid<<".\n";
    if(detailnr>1){
      cout << "Prescribed velocity:" << u_p() << "\n";
      cout << "Prescribed fluid temperature:" << T0 << "\n";      
      // cout << "Prescribed solid temperature:" << prescribeTs.getVals()() << "\n";      
    }
  } // show

} // namespace tube

//////////////////////////////////////////////////////////////////////
