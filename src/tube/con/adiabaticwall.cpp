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
  void AdiabaticWall::show(us i) const {
    TRACE(5,"AdiabaticWall::show()");
    if(isInit()){
      string side;
      if(position==pos::left)
        side="left";
      else
        side="right";
      cout << "AdiabaticWall boundary condition set at << "<<side <<" side of segment .\n";
    }
    else{
      cout << "Show called but init not yet done!\n";
    }
  }

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
    us Ns=gc->Ns();
    error.subvec(0,Ns-1)=Uiszero.error();
    const TubeBcVertex* vertex;
    if(position==pos::left){
      vertex=&t->leftVertex();
    }
    else{
      vertex=&t->rightVertex();
    }    
    error.subvec(Ns,2*Ns-1)=vertex->extrapolateQuant(physquant::heatFlow);
    
    return error;
  }
  void AdiabaticWall::jac(Jacobian& jac) const {
    TRACE(4,"AdiabaticWall::jac()");
    jac+=Uiszero.jac();
    us Ns=gc->Ns();
    JacRow heatFlowjac(firsteqnr+Ns,2);
    const TubeBcVertex* vertex;
    if(position==pos::left){
      vertex=&t->leftVertex();
    }
    else{
      vertex=&t->rightVertex();
    }    
    heatFlowjac+=vertex->dExtrapolateQuant(physquant::heatFlow);
    jac+=heatFlowjac;
  }
} // namespace tube
