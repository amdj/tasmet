#include "adiabaticwall.h"
#include "jacobian.h"
#include "tasystem.h"
#include "bccell.h"
#include "tube.h"
#include "isentropictube.h"
#define Ns (gc->Ns())

namespace tube{

  using tasystem::var;
  using tasystem::TaSystem;
  using tasystem::Globalconf;
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;
  AdiabaticWall::AdiabaticWall(us segnr,Pos position,bool arbitrateMass):
      TubeBc(segnr,position),
      arbitrateMass(arbitrateMass)
  {}

  AdiabaticWall::AdiabaticWall(const AdiabaticWall& o,const TaSystem& sys):
    TubeBc(o,sys),
    arbitrateMass(o.arbitrateMass)
  {
    TRACE(15,"AdiabaticWall::init()");
    if(dynamic_cast<const IsentropicTube*>(t)){
      TRACE(40,"Tube is isentropic");
    }
    setInit(true);
  }
  void AdiabaticWall::updateNf(){
    TRACE(15,"AdiabaticWall::updateNf()");
    massflowzero.updateNf();
    enthalpyflowzero.updateNf(); 
  }
  void AdiabaticWall::show(us i) const {
    TRACE(5,"AdiabaticWall::show()");
    checkInit();
    cout << "AdiabaticWall boundary condition set at "<< posWord(pos) <<" side of segment "<<segnr<<".\n";
  }
  int AdiabaticWall::arbitrateMassEq() const {
    TRACE(15,"AdiabaticWall::arbitrateMassEq()");
    if(arbitrateMass)
      return firsteqnr;
    else
      return -1;
  }
  segment::Connector* AdiabaticWall::copy(const tasystem::TaSystem& sys) const {
    return new AdiabaticWall(*this,sys);
  }
  us AdiabaticWall::getNEqs() const {
    TRACE(15,"AdiabaticWall::getNEqs()");
    return 3*Ns;
  }
  void AdiabaticWall::setEqNrs(us firsteqnr){
    TRACE(2,"AdiabaticWall::setEqNrs()");
    this->firsteqnr=firsteqnr;
    massflowzero.set(firsteqnr,t->bcCell(pos).mbc());
    enthalpyflowzero.set(firsteqnr+Ns,t->bcCell(pos).mHbc());
  }
  vd AdiabaticWall::error() const {
    TRACE(4,"AdiabaticWall::error()");
    const BcCell& cell=t->bcCell(pos);
    vd error(getNEqs());
    error.subvec(0,Ns-1)=massflowzero.error();
    error.subvec(Ns,2*Ns-1)=enthalpyflowzero.error();
    error.subvec(2*Ns,3*Ns-1)=cell.extrapolateQuant(Varnr::Q);
    return error;
  }
  void AdiabaticWall::jac(Jacobian& jac) const {
    TRACE(4,"AdiabaticWall::jac()");
    const BcCell& cell=t->bcCell(pos);
    jac+=massflowzero.jac();
    jac+=enthalpyflowzero.jac();
    JacRow heatflowjac(firsteqnr+2*Ns,2);
    heatflowjac+=cell.dExtrapolateQuant(Varnr::Q);
    jac+=heatflowjac;
  }
} // namespace tube
