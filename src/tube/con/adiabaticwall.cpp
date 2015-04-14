#include "adiabaticwall.h"
#include "jacobian.h"
#include "tasystem.h"
#include "bccell.h"
#include "tube.h"
#include "isentropictube.h"
#define Ns (gc->Ns())

namespace tube{

  using variable::var;
  using tasystem::TaSystem;
  using tasystem::Globalconf;
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  AdiabaticWall::AdiabaticWall(const AdiabaticWall& o,const TaSystem& sys):
    TubeBc(o,sys)
  {
    TRACE(15,"AdiabaticWall::init()");
    if(dynamic_cast<const IsentropicTube*>(t)){
      isentropic=true;
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
    const char* side;
    if(pos==Pos::left)
      side="left";
    else
      side="right";
    cout << "AdiabaticWall boundary condition set at "<<side <<" side of segment "<<segnr<<".\n";
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
    error.subvec(2*Ns,3*Ns-1)=cell.extrapolateQuant(Physquant::HeatFlow);
    return error;
  }
  void AdiabaticWall::jac(Jacobian& jac) const {
    TRACE(4,"AdiabaticWall::jac()");
    const BcCell& cell=t->bcCell(pos);
    jac+=massflowzero.jac();
    jac+=enthalpyflowzero.jac();
    JacRow heatflowjac(firsteqnr+2*Ns,2);
    heatflowjac+=cell.dExtrapolateQuant(Physquant::HeatFlow);
    jac+=heatflowjac;
  }
} // namespace tube
