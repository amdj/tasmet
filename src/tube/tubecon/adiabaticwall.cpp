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
  AdiabaticWall::AdiabaticWall(const string& segid,Pos position,bool arbitrateMass):
      TubeBc(segid,position),
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
  }
  void AdiabaticWall::updateNf(){
    TRACE(15,"AdiabaticWall::updateNf()");
    massflowzero.updateNf();
    enthalpyflowzero.updateNf(); 
  }
  void AdiabaticWall::show(us i) const {
    TRACE(5,"AdiabaticWall::show()");
    cout << "AdiabaticWall boundary condition set at "<< posWord(pos) <<" side of segment "<<segid<<".\n";
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
    if(t->hasSolid())
      return 4*Ns;
    else
      return 3*Ns;
  }
  void AdiabaticWall::setEqNrs(us firsteqnr){
    TRACE(2,"AdiabaticWall::setEqNrs()");
    TubeBc::setEqNrs(firsteqnr);
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
    if(t->hasSolid())
      error.subvec(3*Ns,4*Ns-1)=cell.extrapolateQuant(Varnr::Qs);
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
    // If solid is present, set this heat flow also to zero
    if(t->hasSolid()){
      JacRow solidheatflowjac(firsteqnr+3*Ns,2);
      solidheatflowjac+=cell.dExtrapolateQuant(Varnr::Qs);
      jac+=solidheatflowjac;
    }
  }
} // namespace tube
