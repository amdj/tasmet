#include "adiabaticwall.h"
#include "jacobian.h"
#include "tasystem.h"
#include "bccell.h"
#include "tube.h"
#include "isentropictube.h"


namespace tube{

  using variable::var;
  using tasystem::TaSystem;
  using tasystem::Globalconf;
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  void AdiabaticWall::init(const TaSystem& sys){
    TRACE(15,"AdiabaticWall::init()");
    TubeBc::init(sys);
    try{
      if(dynamic_cast<const IsentropicTube*>(t))
        isentropic=true;
      TRACE(40,"Tube is isentropic");
    }
    catch(std::bad_cast&){}
  }
  void AdiabaticWall::updateNf(){
    TRACE(15,"AdiabaticWall::updateNf()");
    massflowzero.updateNf();
    
  }
  void AdiabaticWall::show(us i) const {
    TRACE(5,"AdiabaticWall::show()");
    checkInit();
    string side;
    if(pos==Pos::left)
      side="left";
    else
      side="right";
    cout << "AdiabaticWall boundary condition set at "<<side <<" side of segment "<<segnr<<".\n";
  }
  us AdiabaticWall::getNEqs() const {
    TRACE(15,"AdiabaticWall::getNEqs()");
    return gc->Ns();
  }
  void AdiabaticWall::setEqNrs(us firsteqnr){
    TRACE(2,"AdiabaticWall::setEqNrs()");
    massflowzero.set(firsteqnr,t->bcCell(pos).mbc());
  }
  vd AdiabaticWall::error() const {
    TRACE(4,"AdiabaticWall::error()");
    return massflowzero.error();
  }
  void AdiabaticWall::jac(Jacobian& jac) const {
    TRACE(4,"AdiabaticWall::jac()");
    jac+=massflowzero.jac();
  }
} // namespace tube
