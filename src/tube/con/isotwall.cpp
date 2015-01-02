#include "isotwall.h"
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

  IsoTWall::IsoTWall(us segnr,pos position,const var& Tbc,const var& Tsbc):
    TubeBc(segnr,position),
    Tbc(Tbc),
    Tsbc(Tsbc)
  {
    TRACE(15,"IsoTWall::IsoTWall()");
  }
  IsoTWall::IsoTWall(const IsoTWall& o):
    TubeBc(o),
    Tbc(o.Tbc),
    Tsbc(o.Tsbc)
  {
    TRACE(15,"IsoTWall::IsoTWall(copy)");
  }

  bool IsoTWall::init(const TaSystem& sys){
    TRACE(15,"IsoTWall::init()");
    if(!TubeBc::init(sys))
      return false;

    Tbc.setGc(*gc);
    Tsbc.setGc(*gc);
    return true;
  }
  void IsoTWall::updateNf(){
    Tbc.updateNf();
    Tsbc.updateNf();
  }
  void IsoTWall::show(us i) const {
    TRACE(5,"IsoTWall::show()");
    if(isInit()){
      string side;
      if(position==pos::left)
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
    if(position==pos::left){
      const TubeVertex& vertex=t->leftVertex();
      Uiszero.set(firsteqnr,vertex.UL(),zero);
      Tbc.set(firsteqnr+Ns,vertex.TL()); // vals are set in constructor
      Tsbc.set(firsteqnr+2*Ns,vertex.TsL()); // idem
    }
    else{
      const TubeVertex& vertex=t->rightVertex();
      Uiszero.set(firsteqnr,vertex.UR(),zero);
      Tbc.set(firsteqnr+Ns,vertex.TR()); // vals are set in constructor
      Tsbc.set(firsteqnr+2*Ns,vertex.TsR()); // idem

    }
  }
  vd IsoTWall::error() const {
    TRACE(15,"IsoTWall::error()");
    vd error(getNEqs());
    us Ns=gc->Ns();
    vd rhoextrapolate;
    const TubeBcVertex* vertex;
    if(position==pos::left){
      vertex=&t->leftVertex();
      rhoextrapolate=vertex->rhoL()();
    }
    else{
      vertex=&t->rightVertex();
      rhoextrapolate=vertex->rhoR()();
    }    
    rhoextrapolate-=vertex->extrapolateQuant(physquant::density);

    error.subvec(0,Ns-1)=Uiszero.error();
    error.subvec(Ns,2*Ns-1)=Tbc.error();
    error.subvec(2*Ns,3*Ns-1)=Tsbc.error();
    error.subvec(3*Ns,4*Ns-1)=rhoextrapolate;

    return error;
  }
  void IsoTWall::jac(Jacobian& jac) const {
    TRACE(15,"IsoTWall::jac()");
    us Ns=gc->Ns();
    const TubeBcVertex* vertex;

    JacRow extrapolaterho(firsteqnr+3*Ns,3);
    extrapolaterho.setRowDof(firsteqnr+3*Ns);
    if(position==pos::left){
      vertex=&t->leftVertex();
      extrapolaterho+=JacCol(vertex->rhoL(),eye<dmat>(gc->Ns(),gc->Ns()));
    }
    else{
      vertex=&t->rightVertex();
      extrapolaterho+=JacCol(vertex->rhoR(),eye<dmat>(gc->Ns(),gc->Ns()));
    }    
    extrapolaterho+=(vertex->dExtrapolateQuant(physquant::density)*=-1);
    // Put them into jacobian
    jac+=Uiszero.jac();
    jac+=Tbc.jac();
    jac+=Tsbc.jac();
    jac+=extrapolaterho;
  }
} // namespace tube

