#define TRACERPLUS 20

#include "tube.h"
#include "pressurebc.h"
#include "tasystem.h"
#include "globalconf.h"
#include "jacobian.h"
#include "tubebcvertex.h"
#include "varnr.h"
#include "state.h"

namespace tube{
  using variable::var;
  using tasystem::TaSystem;
  using tasystem::Globalconf;
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;
  
  PressureBc::PressureBc(const var& pres,const var& temp,const var& stemp,us segnr,pos position):
    TubeBc(segnr,position),
    pbc(pres),
    Tbc(temp),
    Tsbc(stemp)
  {
    TRACE(8,"PressureBc full constructor");
  }
  PressureBc::PressureBc(const var& pres,const var& temp,us segnr,pos position):
    PressureBc(pres,temp,coldtemp(pres),segnr,position)
  {}
  PressureBc::PressureBc(const var& pres,us segnr,pos position):
    PressureBc(pres,adiabatictemp(pres),segnr,position)
  {}
  PressureBc::PressureBc(const PressureBc& other):
    TubeBc(other),
    pbc(other.pbc),
    Tbc(other.Tbc),
    Tsbc(other.Tsbc)
  {
    TRACE(8,"PressureBc copy constructor");
  }
  var PressureBc::adiabatictemp(const var& pres){
    TRACE(10,"PressureBc::adiabatictemp()");
    const Globalconf& gc=pres.gc();
    d T0=gc.T0;
    d gamma=gc.gas.gamma(T0);
    vd p0(gc.Ns(),fillwith::ones); p0*=gc.p0;
    vd Tbct=T0*pow((p0+pres.tdata())/p0,(gamma-1.0)/gamma);		// Adiabatic compression/expansion
    var res(pres.gc());
    res.settdata(Tbct);
    return res;
  }
  var PressureBc::coldtemp(const var& pres){
    vd adata(pres.gc().Ns(),fillwith::zeros);
    adata(0)=pres.gc().T0;
    var res(pres.gc());
    res.set(adata);
    return res;
  }
  void PressureBc::updateNf(){
    pbc.updateNf();
    Tbc.updateNf();
    Tsbc.updateNf();
  }

  void PressureBc::init(const TaSystem& sys)
  {
    TRACE(8,"PressureBc::init()");
    TubeBc::init(sys);

    // Decouple from old globalconf pointer
    pbc.setGc(sys.gc);
    Tbc.setGc(sys.gc); 
    Tsbc.setGc(sys.gc);
    TRACE(15,"PressureBc::init()");
  }
  void PressureBc::setEqNrs(us firsteqnr){
    TRACE(2,"Pressure::setEqNrs()");
    this->firsteqnr=firsteqnr;
    us Ns=gc->Ns();
    if(position==pos::left){
      const TubeVertex& vertex=t->leftVertex();
      prescribep.set(firsteqnr,vertex.pL(),pbc());
      prescribeT.set(firsteqnr+Ns,vertex.TL(),Tbc());
      prescribeTs.set(firsteqnr+2*Ns,vertex.TsL(),Tsbc());
    }
    else{
      const TubeVertex& vertex=t->rightVertex();
      prescribep.set(firsteqnr,vertex.pR(),pbc());
      prescribeT.set(firsteqnr+Ns,vertex.TR(),Tbc());
      prescribeTs.set(firsteqnr+2*Ns,vertex.TsR(),Tsbc());
    }
  }
  vd PressureBc::error() const {
    TRACE(15,"PressureBc::error()");
    vd error(getNEqs());
    us Ns=gc->Ns();
    const dmat& fDFT=gc->fDFT;
    vd massflowv;
    const TubeBcVertex* vertex;
    if(position==pos::left){
      vertex=&t->leftVertex();
      massflowv=fDFT*(vertex->rhoL().tdata()%vertex->UL().tdata());
    }
    else{
      vertex=&t->rightVertex();
      massflowv=fDFT*(vertex->rhoR().tdata()%vertex->UR().tdata());
    }

    error.subvec(0,gc->Ns()-1)=prescribep.error();
    error.subvec(Ns,2*Ns-1)=prescribeT.error();
    error.subvec(2*Ns,3*Ns-1)=prescribeTs.error();
    error.subvec(3*Ns,4*Ns-1)=massflowv-vertex->extrapolateQuant(physquant::massFlow);

    VARTRACE(15,error)
    return error;
  }
  void PressureBc::jac(Jacobian& jac) const{
    TRACE(8,"PressureBc::jac()");
    us Ns=gc->Ns();
    JacRow massfloweq(firsteqnr+3*Ns,6);

    const dmat& fDFT=gc->fDFT;
    const dmat& iDFT=gc->iDFT;

    const TubeBcVertex* vertex;
    if(position==pos::left){
      vertex=&t->leftVertex();
      massfloweq+=JacCol(vertex->rhoL(),fDFT*vertex->UL().diagt()*iDFT);
      massfloweq+=JacCol(vertex->UL(),fDFT*vertex->rhoL().diagt()*iDFT);
    }
    else{
      vertex=&t->rightVertex();
      massfloweq+=JacCol(vertex->rhoR(),fDFT*vertex->UR().diagt()*iDFT);
      massfloweq+=JacCol(vertex->UR(),fDFT*vertex->rhoR().diagt()*iDFT);
    }
    massfloweq+=(vertex->dExtrapolateQuant(physquant::massFlow)*=-1);

    jac+=prescribep.jac();
    jac+=prescribeT.jac();
    jac+=prescribeTs.jac();
    jac+=massfloweq;

  }

  void PressureBc::show(us i) const {
    TRACE(5,"PressureBc::show()");
    if(isInit()){
      string side;
      if(position==pos::left)
        side="left";
      else
        side="right";
      cout << "PressureBc boundary condition set at << "<<side <<" side of segment .\n";
      cout << "Prescribed pressure:" << pbc() << "\n";
    }
    else{
      cout << "Show called but init not yet done!\n";
    }
  }



} // namespace tube



