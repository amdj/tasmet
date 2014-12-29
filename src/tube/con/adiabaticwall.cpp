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
    zero=vd(sys.gc.Ns(),fillwith::zeros);
    return true;
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
      // d xR=vertex.weightFactors().xR;
      // d xRR=vertex.right()->weightFactors().xR;
      // dpdxiszero.set(firsteqnr+Ns,vertex.pL(),vertex.pR(),vertex.right()->pR(),0,xR,xRR,zero);
      d xR=vertex.weightFactors().vx;
      d xRR=vertex.right()->weightFactors().vx;
      dpdxiszero.set(firsteqnr+Ns,vertex.rhoL(),vertex.rho(),vertex.right()->rho(),0,xR,xRR,zero);
    }
    else{
      const TubeVertex& vertex=t->rightVertex();
      Uiszero.set(firsteqnr,vertex.UR(),zero);
      // d xi=vertex.weightFactors().xR;
      // d xj=vertex.weightFactors().xL;
      // d xk=vertex.left()->weightFactors().xL;
      // dpdxiszero.set(firsteqnr+Ns,vertex.pR(),vertex.pL(),vertex.left()->pL(),xi,xj,xk,zero);
      d xi=vertex.weightFactors().xR;
      d xj=vertex.weightFactors().vx;
      d xk=vertex.left()->weightFactors().vx;
      dpdxiszero.set(firsteqnr+Ns,vertex.rhoR(),vertex.rho(),vertex.left()->rho(),xi,xj,xk,zero);
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
    error.subvec(Ns,2*Ns-1)=dpdxiszero.error();
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

    JacRow dpdx(firsteqnr+3*Ns,2);
    dpdx+=JacCol(vertex->pL(),eye(Ns,Ns));
    dpdx+=JacCol(vertex->pR(),-eye(Ns,Ns));


    // Put them into jacobian
    jac+=Uiszero.jac();
    jac+=dpdxiszero.jac();
    jac+=heatFlowjac;
    jac+=solidheatFlowjac;
    // JacRow massfloweq(firsteqnr+3*Ns,4);    
    // massfloweq+=vertex->dExtrapolateQuant(physquant::massFlow);
    // jac+=massfloweq;
  }
} // namespace tube
