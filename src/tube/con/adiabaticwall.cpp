#include "adiabaticwall.h"
#include "jacobian.h"
#include "tasystem.h"
#include "tubebccell.h"
#include "tube.h"
#include "weightfactors.h"

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
    Uiszero.setGc(*gc); 
    drhodxiszero.setGc(*gc);
    setInit(true);
  }
  void AdiabaticWall::updateNf(){
    TRACE(15,"AdiabaticWall::updateNf()");
    Uiszero.updateNf();
    drhodxiszero.updateNf();
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

  void AdiabaticWall::setEqNrs(us firsteqnr){
    TRACE(2,"AdiabaticWall::setEqNrs()");
    this->firsteqnr=firsteqnr;
    us Ns=gc->Ns();

    var zero(*gc,0);             // zero at all the time

    if(pos==Pos::left){
      const Cell& cell=t->leftCell();
      Uiszero.set(firsteqnr,cell.UL(),zero);
      // d xR=cell.weightFactors().xR;
      // d xRR=cell.right()->weightFactors().xR;
      // drhodxiszero.set(firsteqnr+Ns,cell.pL(),cell.pR(),cell.right()->pR(),0,xR,xRR,zero);
      d xR=cell.weightFactors().vx;
      d xRR=cell.right()->weightFactors().vx;
      drhodxiszero.set(firsteqnr+Ns,cell.rhoL(),cell.rho(),cell.right()->rho(),0,xR,xRR,zero);
    }
    else{
      const Cell& cell=t->rightCell();
      Uiszero.set(firsteqnr,cell.UR(),zero);
      // d xi=cell.weightFactors().xR;
      // d xj=cell.weightFactors().xL;
      // d xk=cell.left()->weightFactors().xL;
      // drhodxiszero.set(firsteqnr+Ns,cell.pR(),cell.pL(),cell.left()->pL(),xi,xj,xk,zero);
      d xi=cell.weightFactors().xR;
      d xj=cell.weightFactors().vx;
      d xk=cell.left()->weightFactors().vx;
      drhodxiszero.set(firsteqnr+Ns,cell.rhoR(),cell.rho(),cell.left()->rho(),xi,xj,xk,zero);
    }
  }
  vd AdiabaticWall::error() const {
    TRACE(4,"AdiabaticWall::error()");
    vd error(getNEqs());
    us Ns=gc->Ns();
    const TubeBcCell* cell;
    if(pos==Pos::left){
      cell=&t->leftCell();
    }
    else{
      cell=&t->rightCell();
    }    
    error.subvec(0,Ns-1)=Uiszero.error();
    error.subvec(Ns,2*Ns-1)=drhodxiszero.error();
    error.subvec(2*Ns,3*Ns-1)=cell->extrapolateQuant(physquant::heatFlow);
    error.subvec(3*Ns,4*Ns-1)=cell->extrapolateQuant(physquant::solidHeatFlow);

    return error;
  }
  void AdiabaticWall::jac(Jacobian& jac) const {
    TRACE(4,"AdiabaticWall::jac()");
    us Ns=gc->Ns();
    const TubeBcCell* cell;
    if(pos==Pos::left){
      cell=&t->leftCell();
    }
    else{
      cell=&t->rightCell();
    }    
    JacRow heatFlowjac=cell->dExtrapolateQuant(physquant::heatFlow);
    heatFlowjac.setRowDof(firsteqnr+2*Ns);
    JacRow solidheatFlowjac=cell->dExtrapolateQuant(physquant::solidHeatFlow);
    solidheatFlowjac.setRowDof(firsteqnr+3*Ns);

    // Put them into jacobian
    jac+=Uiszero.jac();
    jac+=drhodxiszero.jac();
    jac+=heatFlowjac;
    jac+=solidheatFlowjac;
  }
} // namespace tube
