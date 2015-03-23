#include "adiabaticwall.h"
#include "jacobian.h"
#include "tasystem.h"
#include "tubebcvertex.h"
#include "tube.h"
#include "weightfactors.h"

namespace tube{

  using variable::var;
  using tasystem::TaSystem;
  using tasystem::Globalconf;
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  bool AdiabaticWall::init(const TaSystem& sys){
    TRACE(15,"AdiabaticWall::init()");
    if(!TubeBc::init(sys))
      return false;
    Uiszero.setGc(*gc); 
    drhodxiszero.setGc(*gc);
    setInit(true);
    return true;
  }
  void AdiabaticWall::updateNf(){
    TRACE(15,"AdiabaticWall::updateNf()");
    Uiszero.updateNf();
    drhodxiszero.updateNf();
  }
  void AdiabaticWall::show(us i) const {
    TRACE(5,"AdiabaticWall::show()");
    if(isInit()){
      string side;
      if(position==pos::left)
        side="left";
      else
        side="right";
      cout << "AdiabaticWall boundary condition set at "<<side <<" side of segment "<<segnr<<".\n";
    }
    else
      WARN("Show called but init not yet done!");
  }

  void AdiabaticWall::setEqNrs(us firsteqnr){
    TRACE(2,"AdiabaticWall::setEqNrs()");
    this->firsteqnr=firsteqnr;
    us Ns=gc->Ns();

    var zero(*gc,0);             // zero at all the time

    if(position==pos::left){
      const TubeVertex& vertex=t->leftVertex();
      Uiszero.set(firsteqnr,vertex.UL(),zero);
      // d xR=vertex.weightFactors().xR;
      // d xRR=vertex.right()->weightFactors().xR;
      // drhodxiszero.set(firsteqnr+Ns,vertex.pL(),vertex.pR(),vertex.right()->pR(),0,xR,xRR,zero);
      d xR=vertex.weightFactors().vx;
      d xRR=vertex.right()->weightFactors().vx;
      drhodxiszero.set(firsteqnr+Ns,vertex.rhoL(),vertex.rho(),vertex.right()->rho(),0,xR,xRR,zero);
    }
    else{
      const TubeVertex& vertex=t->rightVertex();
      Uiszero.set(firsteqnr,vertex.UR(),zero);
      // d xi=vertex.weightFactors().xR;
      // d xj=vertex.weightFactors().xL;
      // d xk=vertex.left()->weightFactors().xL;
      // drhodxiszero.set(firsteqnr+Ns,vertex.pR(),vertex.pL(),vertex.left()->pL(),xi,xj,xk,zero);
      d xi=vertex.weightFactors().xR;
      d xj=vertex.weightFactors().vx;
      d xk=vertex.left()->weightFactors().vx;
      drhodxiszero.set(firsteqnr+Ns,vertex.rhoR(),vertex.rho(),vertex.left()->rho(),xi,xj,xk,zero);
    }
  }
  vd AdiabaticWall::error() const {
    TRACE(4,"AdiabaticWall::error()");
    vd error(getNEqs());
    us Ns=gc->Ns();
    const TubeBcVertex* vertex;
    if(position==pos::left){
      vertex=&t->leftVertex();
    }
    else{
      vertex=&t->rightVertex();
    }    
    error.subvec(0,Ns-1)=Uiszero.error();
    error.subvec(Ns,2*Ns-1)=drhodxiszero.error();
    error.subvec(2*Ns,3*Ns-1)=vertex->extrapolateQuant(physquant::heatFlow);
    error.subvec(3*Ns,4*Ns-1)=vertex->extrapolateQuant(physquant::solidHeatFlow);

    return error;
  }
  void AdiabaticWall::jac(Jacobian& jac) const {
    TRACE(4,"AdiabaticWall::jac()");
    us Ns=gc->Ns();
    const TubeBcVertex* vertex;
    if(position==pos::left){
      vertex=&t->leftVertex();
    }
    else{
      vertex=&t->rightVertex();
    }    
    JacRow heatFlowjac=vertex->dExtrapolateQuant(physquant::heatFlow);
    heatFlowjac.setRowDof(firsteqnr+2*Ns);
    JacRow solidheatFlowjac=vertex->dExtrapolateQuant(physquant::solidHeatFlow);
    solidheatFlowjac.setRowDof(firsteqnr+3*Ns);

    // Put them into jacobian
    jac+=Uiszero.jac();
    jac+=drhodxiszero.jac();
    jac+=heatFlowjac;
    jac+=solidheatFlowjac;
  }
} // namespace tube
