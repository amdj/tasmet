#include "adiabaticwall.h"
#include "jacobian.h"
#include "tasystem.h"
#include "tubebcvertex.h"
#include "tube.h"

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
    zero=vd(sys.gc.Ns(),fillwith::zeros);
  }
  void AdiabaticWall::updateNf(){}
  void AdiabaticWall::setEqNrs(us firsteqnr){
    TRACE(2,"AdiabaticWall::setEqNrs()");
    this->firsteqnr=firsteqnr;
    us Ns=gc->Ns();
    if(position==pos::left){
      const TubeVertex& vertex=t->leftVertex();
      Uiszero.set(firsteqnr,vertex.UL(),zero);
    }
    else{
      const TubeVertex& vertex=t->rightVertex();
      Uiszero.set(firsteqnr,vertex.UR(),zero);
    }
  }
  vd AdiabaticWall::error() const {
    TRACE(4,"AdiabaticWall::error()");
    vd error(getNEqs());
    error.subvec(0,gc->Ns()-1)=Uiszero.error();
  }
  void AdiabaticWall::jac(Jacobian& jac) const {
    TRACE(4,"AdiabaticWall::jac()");
    jac+=Uiszero.jac();
    
  }
} // namespace tube
