#define TRACERPLUS 20

#include "tube.h"
#include "pressurebc.h"
#include "tasystem.h"
#include "globalconf.h"
#include "jacobian.h"
#include "tubebcvertex.h"
#include "constants.h"
#include "state.h"

namespace tube{
  using variable::var;
  using tasystem::TaSystem;
  using tasystem::Globalconf;
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  var coldtemp(const var& pres){
    return var(pres.gc(),pres.gc().T0());
  }
  
  PressureBc::PressureBc(const var& pres,const var& temp,const var& stemp,us segnr,Pos position):
    TubeBc(segnr,position),
    prescribep(pres),
    prescribeT(temp),
    prescribeTs(stemp)
  {
    TRACE(8,"PressureBc full constructor");
  }
  PressureBc::PressureBc(const var& pres,const var& temp,us segnr,Pos position):
    PressureBc(pres,temp,coldtemp(pres),segnr,position)
  {}
  PressureBc::PressureBc(const var& pres,us segnr,Pos position):
    PressureBc(pres,adiabatictemp(pres),segnr,position)
  {}
  PressureBc::PressureBc(const PressureBc& other):
    TubeBc(other),
    prescribep(other.prescribep),
    prescribeT(other.prescribeT),
    prescribeTs(other.prescribeTs)
  {
    TRACE(8,"PressureBc copy constructor");
  }
  var PressureBc::adiabatictemp(const var& pres){
    TRACE(10,"PressureBc::adiabatictemp()");
    const Globalconf& gc=pres.gc();
    d T0=gc.T0();
    d gamma=gc.gas().gamma(T0);
    vd p0(gc.Ns(),fillwith::ones); p0*=gc.p0();
    vd Tbct=T0*pow((p0+pres.tdata())/p0,(gamma-1.0)/gamma);		// Adiabatic compression/expansion
    var res(pres.gc());
    res.settdata(Tbct);
    return res;
  }
  void PressureBc::updateNf(){
    prescribep.updateNf();
    prescribeT.updateNf();
    prescribeTs.updateNf();
  }
  void PressureBc::init(const TaSystem& sys)
  {
    TRACE(8,"PressureBc::init()");
    TubeBc::init(sys);

    // Decouple from old globalconf pointer
    prescribep.setGc(*gc);
    prescribeT.setGc(*gc);
    prescribeTs.setGc(*gc); 
    setInit(true);
  }
  void PressureBc::setEqNrs(us firsteqnr){
    TRACE(2,"Pressure::setEqNrs()");
    this->firsteqnr=firsteqnr;
    us Ns=gc->Ns();
    if(pos==Pos::left){
      const TubeVertex& vertex=t->leftVertex();
      prescribep.set(firsteqnr,vertex.pL());
      prescribeT.set(firsteqnr+Ns,vertex.TL());
      prescribeTs.set(firsteqnr+2*Ns,vertex.TsL());
    }
    else{
      const TubeVertex& vertex=t->rightVertex();
      prescribep.set(firsteqnr,vertex.pR());
      prescribeT.set(firsteqnr+Ns,vertex.TR());
      prescribeTs.set(firsteqnr+2*Ns,vertex.TsR());
    }
  }
  vd PressureBc::error() const {
    TRACE(15,"PressureBc::error()");
    vd error(getNEqs());
    us Ns=gc->Ns();
    const dmat& fDFT=gc->fDFT;
    vd massflowv;
    const TubeBcVertex* vertex;
    if(pos==Pos::left){
      vertex=&t->leftVertex();
      massflowv=fDFT*(vertex->rhoL().tdata()%vertex->UL().tdata());
    }
    else{
      vertex=&t->rightVertex();
      massflowv=fDFT*(vertex->rhoR().tdata()%vertex->UR().tdata());
    }

    error.subvec(0,gc->Ns()-1)=prescribep.error();
    error.subvec(1*Ns,2*Ns-1)=prescribeT.error();
    error.subvec(2*Ns,3*Ns-1)=prescribeTs.error();
    error.subvec(3*Ns,4*Ns-1)=massflowv-vertex->extrapolateQuant(physquant::massFlow);

    // VARTRACE(15,error)
    return error;
  }
  void PressureBc::jac(Jacobian& jac) const{
    TRACE(8,"PressureBc::jac()");
    us Ns=gc->Ns();
    JacRow massfloweq(firsteqnr+3*Ns,6);

    const dmat& fDFT=gc->fDFT;
    const dmat& iDFT=gc->iDFT;

    const TubeBcVertex* vertex;
    if(pos==Pos::left){
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

  void PressureBc::show(us detailnr) const {
    TRACE(5,"PressureBc::show()");
    if(isInit()){
      string side;
      if(pos==Pos::left)
        side="left";
      else
        side="right";
      cout << "PressureBc boundary condition set at "<<side <<" side of segment "<< segnr<<".\n";
      if(detailnr>1){
        cout << "Prescribed pressure:" << prescribep.getVals()() << "\n";
        cout << "Prescribed fluid temperature:" << prescribeT.getVals()() << "\n";      
        cout << "Prescribed solid temperature:" << prescribeTs.getVals()() << "\n";      
      }
    }
    else{
      WARN("Show called but init not yet done!");
    } // isInit()
  } // show

} // namespace tube



