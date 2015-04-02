#include "tube.h"
#include "pressurebc.h"
#include "tasystem.h"
#include "globalconf.h"
#include "jacobian.h"
#include "bccell.h"
#include "constants.h"
#include "state.h"

#define fDFT (gc->fDFT)
#define iDFT (gc->iDFT)
#define DDTfd (gc->DDTfd)

#define Ns (gc->Ns())

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
  var adiabatictemp(const var& pres){
    TRACE(10,"adiabatictemp()");
    const Globalconf* gc=&pres.gc();
    d T0=gc->T0();
    d gamma=gc->gas().gamma(T0);
    vd p0(Ns,fillwith::ones); p0*=gc->p0();
    vd Tbct=T0*pow((p0+pres.tdata())/p0,(gamma-1.0)/gamma);		// Adiabatic compression/expansion
    var res(pres.gc());
    res.settdata(Tbct);
    return res;
  } 
  PressureBc::PressureBc(const var& pres,const var& temp,const var& stemp,us segnr,Pos position):
    TubeBc(segnr,position),
    p_prescribed(pres),
    prescribeT(temp)
  {
    
    TRACE(8,"PressureBc full constructor");
    VARTRACE(30,temp());
  }
  PressureBc::PressureBc(const var& pres,const var& temp,us segnr,Pos position):
    PressureBc(pres,temp,coldtemp(pres),segnr,position)
  {}
  PressureBc::PressureBc(const var& pres,us segnr,Pos position):
    PressureBc(pres,adiabatictemp(pres),segnr,position)
  {}
  PressureBc::PressureBc(const PressureBc& other):
    TubeBc(other),
    p_prescribed(other.p_prescribed),
    prescribeT(other.prescribeT)
  {
    TRACE(8,"PressureBc copy constructor");
  }
 
  void PressureBc::updateNf(){
    p_prescribed.updateNf();
    prescribeT.updateNf();
  }
  void PressureBc::init(const TaSystem& sys)
  {
    TRACE(8,"PressureBc::init()");
    TubeBc::init(sys);
    assert(gc);

    // Decouple from old globalconf pointer
    p_prescribed.setGc(*gc);
    prescribeT.setGc(*gc);
    setInit(true);
  }
  void PressureBc::setEqNrs(us firsteqnr){
    TRACE(2,"Pressure::setEqNrs()");
    const BcCell& cell=t->bcCell(pos);
    this->firsteqnr=firsteqnr;
    prescribeT.set(firsteqnr+Ns,cell.Tbc());
  }
  vd PressureBc::error() const {
    TRACE(15,"PressureBc::error()");
    vd error(getNEqs(),fillwith::zeros);
    const BcCell& cell=t->bcCell(pos);

    vd errorM(Ns,fillwith::zeros);
    if(pos==Pos::left)    {
      d Wddt=cell.vx;
      errorM+=Wddt*DDTfd*cell.mbc()();
      errorM+=cell.vSf*cell.p()();
      errorM-=cell.SfL*p_prescribed();
      errorM+=cell.mu()();       // 
      errorM-=cell.extrapolateQuant(Physquant::momentumFlow);
      errorM+=Wddt*(t->getDragResistance().drag(cell));
    }
    else{
      d Wddt=cell.xR-cell.vx;
      errorM+=Wddt*DDTfd*cell.mbc()();
      errorM-=cell.vSf*cell.p()();
      errorM+=cell.SfR*p_prescribed();
      errorM-=cell.mu()();
      errorM+=cell.extrapolateQuant(Physquant::momentumFlow);      
      errorM+=Wddt*(t->getDragResistance().drag(cell));
    }      
    error.subvec(0,Ns-1)=errorM;
    error.subvec(Ns,2*Ns-1)=prescribeT.error();
    return error;
  }
  void PressureBc::jac(Jacobian& jac) const{
    TRACE(8,"PressureBc::jac()");

    const BcCell& cell=t->bcCell(pos);
    if(pos==Pos::left)    {
      d Wddt=cell.vx;
      VARTRACE(25,Wddt);
      JacRow jacr(firsteqnr,4);
      jacr+=JacCol(cell.mbc(),Wddt*DDTfd);
      jacr+=JacCol(cell.p(),cell.vSf*eye(Ns,Ns));
      jacr+=JacCol(cell.mu(),eye(Ns,Ns));
      jacr+=-cell.dExtrapolateQuant(Physquant::momentumFlow);
      jacr+=JacCol(cell.mL(),Wddt*(t->getDragResistance().dm(cell)));
      jac+=jacr;
    }
    else{
      d Wddt=cell.xR-cell.vx;
      VARTRACE(25,Wddt);
      JacRow jacr(firsteqnr,4);
      jacr+=JacCol(cell.mbc(),Wddt*DDTfd);
      jacr+=JacCol(cell.p(),-cell.vSf*eye(Ns,Ns));
      jacr+=JacCol(cell.mu(),-eye(Ns,Ns));
      jacr+=cell.dExtrapolateQuant(Physquant::momentumFlow);
      jacr+=JacCol(cell.mR(),Wddt*(t->getDragResistance().dm(cell)));
      jac+=jacr;
    }
    jac+=prescribeT.jac();
  }
  void PressureBc::show(us detailnr) const {
    TRACE(5,"PressureBc::show()");
    checkInit();
    string side;
    if(pos==Pos::left)
      side="left";
    else
      side="right";
    cout << "PressureBc boundary condition set at "<<side <<" side of segment "<< segnr<<".\n";
    if(detailnr>1){
      cout << "Prescribed pressure:" << p_prescribed() << "\n";
      cout << "Prescribed fluid temperature:" << prescribeT.getVals()() << "\n";      
      // cout << "Prescribed solid temperature:" << prescribeTs.getVals()() << "\n";      
    }
  } // show

} // namespace tube



