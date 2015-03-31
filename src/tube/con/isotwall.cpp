#include "isotwall.h"
#include "jacobian.h"
#include "tasystem.h"
#include "tubebccell.h"
#include "tube.h"
#include "weightfactors.h"

#define iDFT (gc->iDFT)
#define fDFT (gc->fDFT)
#define DDTfd (gc->DDTfd)

namespace tube{

  using variable::var;
  using tasystem::TaSystem;
  using tasystem::Globalconf;
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  IsoTWall::IsoTWall(us segnr,Pos position,const var& Tbc):
    TubeBc(segnr,position),
    Tbc(Tbc),
    Tsbc(Tbc)
  {
    TRACE(15,"IsoTWall::IsoTWall()");
  }
  IsoTWall::IsoTWall(const IsoTWall& o):
    TubeBc(o),
    Uiszero(o.Uiszero),
    Tbc(o.Tbc),
    Tsbc(o.Tsbc)
  {
    TRACE(15,"IsoTWall::IsoTWall(copy)");
  }

  void IsoTWall::init(const TaSystem& sys){
    TRACE(15,"IsoTWall::init()");
    TubeBc::init(sys);
    Uiszero.setGc(*gc);
    Tbc.setGc(*gc);
    Tsbc.setGc(*gc);
  }
  void IsoTWall::updateNf(){
    Tbc.updateNf();
    Tsbc.updateNf();
  }
  void IsoTWall::show(us i) const {
    TRACE(5,"IsoTWall::show()");
    if(isInit()){
      string side;
      if(pos==Pos::left)
        side="left";
      else
        side="right";
      cout << "IsoTWall boundary condition set at << "<<side <<" side of segment .\n";
    }
    else{
      cout << "Show called but init not yet done!\n";
    }
  }
  void IsoTWall::setEqNrs(us firsteqnr){
    TRACE(2,"IsoTWall::setEqNrs()");
    this->firsteqnr=firsteqnr;
    var zero(*gc,0);             // zero at all the time

    us Ns=gc->Ns();
    if(pos==Pos::left){
      const Cell& cell=t->leftCell();
      Uiszero.set(firsteqnr,cell.UL(),zero);
      Tbc.set(firsteqnr+Ns,cell.TL()); // vals are set in constructor
      Tsbc.set(firsteqnr+2*Ns,cell.TsL()); // idem

    }
    else{
      const Cell& cell=t->rightCell();
      Uiszero.set(firsteqnr,cell.UR(),zero);
      Tbc.set(firsteqnr+Ns,cell.TR()); // vals are set in constructor
      Tsbc.set(firsteqnr+2*Ns,cell.TsR()); // idem

    }
  }
  vd IsoTWall::error() const {
    TRACE(15,"IsoTWall::error()");
    vd error(getNEqs());
    us Ns=gc->Ns();

    const TubeBcCell* cell;
    
    const var *pb;
    if(pos==Pos::left){
      cell=&t->leftCell();
      pb=&cell->pL();
    }
    else{
      cell=&t->rightCell();
      pb=&cell->pR();
    }    
    // Compute error on state eq derivative
    vd stateer=(*pb)(); stateer(0)+=gc->p0();
    stateer+=-cell->extrapolateQuant(physquant::rhoRT);
    VARTRACE(30,stateer);
    // Add all individual error parts to vector
    error.subvec(0,Ns-1)=Uiszero.error();
    error.subvec(Ns,2*Ns-1)=Tbc.error();
    error.subvec(2*Ns,3*Ns-1)=Tsbc.error();
    // error.subvec(3*Ns,4*Ns-1)=cell->extrapolateQuant(physquant::massFlow);
    error.subvec(3*Ns,4*Ns-1)=stateer;
    return error;
  }
  void IsoTWall::jac(Jacobian& jac) const {
    TRACE(15,"IsoTWall::jac()");
    us Ns=gc->Ns();

    const TubeBcCell* cell;
    const variable::var* pb;

    if(pos==Pos::left){
      cell=&t->leftCell();
      pb=&cell->pL();
    }
    else{
      cell=&t->rightCell();
      pb=&cell->pR();
    }

    // Derivative extrapolated state equation
    JacRow stateeq(firsteqnr+3*Ns,6);
    stateeq+=JacCol(*pb,eye<dmat>(gc->Ns(),gc->Ns()));
    stateeq+=(cell->dExtrapolateQuant(physquant::rhoRT)*=-1);


    // Put them into jacobian
    jac+=Uiszero.jac();
    jac+=Tbc.jac();
    jac+=Tsbc.jac();
    // jac+=massfloweq;
    jac+=stateeq;
  }
} // namespace tube

