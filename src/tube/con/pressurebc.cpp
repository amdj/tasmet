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
#define eye (eye(Ns,Ns))

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
  namespace 
  {
    inline d cp(const Cell& c) {
      return c.gc->gas().cp(c.gc->T0());
    }
  } // namespace 
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
      #ifndef NODRAG
      errorM+=Wddt*(t->getDragResistance().drag(cell));
      #endif
    }
    else{
      d Wddt=cell.xR-cell.vx;
      errorM+=Wddt*DDTfd*cell.mbc()();
      errorM-=cell.vSf*cell.p()();
      errorM+=cell.SfR*p_prescribed();
      errorM-=cell.mu()();
      errorM+=cell.extrapolateQuant(Physquant::momentumFlow);      
      #ifndef NODRAG
      errorM+=Wddt*(t->getDragResistance().drag(cell));
      #endif
    }      
    error.subvec(0,Ns-1)=errorM;
    error.subvec(Ns,2*Ns-1)=prescribeT.error();

    error.subvec(2*Ns,3*Ns-1)=-cell.mHbc()()
      +fDFT*(cp(cell)*cell.mbc().tdata()%cell.Tbc().tdata());
    // VARTRACE(30,cell.Tbc()());
    // VARTRACE(30,cp(cell));
    // error.subvec(2*Ns,3*Ns-1)=-cell.mHbc()();
    // error.subvec(2*Ns,3*Ns-1)+=cell.extrapolateQuant(Physquant::enthalpyFlow);
    VARTRACE(20,error);
    return error;
  }
  void PressureBc::jac(Jacobian& jac) const{
    TRACE(8,"PressureBc::jac()");

    const BcCell& cell=t->bcCell(pos);
    // Prescribed temperature Jacobian part
    jac+=prescribeT.jac();
    // Momentum equation Jacobian
    if(pos==Pos::left)    {
      d Wddt=cell.vx;
      VARTRACE(25,Wddt);
      JacRow jacr(firsteqnr,4);
      jacr+=JacCol(cell.mbc(),Wddt*DDTfd);
      jacr+=JacCol(cell.p(),cell.vSf*eye);
      jacr+=JacCol(cell.mu(),eye);
      jacr+=-cell.dExtrapolateQuant(Physquant::momentumFlow);
      #ifndef NODRAG
      jacr+=JacCol(cell.mL(),Wddt*(t->getDragResistance().dm(cell)));
      #endif  // NODRAG
      jac+=jacr;
    }
    else{
      d Wddt=cell.xR-cell.vx;
      VARTRACE(25,Wddt);
      JacRow jacr(firsteqnr,4);
      jacr+=JacCol(cell.mbc(),Wddt*DDTfd);
      jacr+=JacCol(cell.p(),-cell.vSf*eye);
      jacr+=JacCol(cell.mu(),-eye);
      jacr+=cell.dExtrapolateQuant(Physquant::momentumFlow);
      #ifndef NODRAG
      jacr+=JacCol(cell.mR(),Wddt*(t->getDragResistance().dm(cell)));
      #endif  // NODRAG
      jac+=jacr;
    }
    // Prescribed enthalpy flow.
    JacRow enthalpy_extrapolated_jac(firsteqnr+2*Ns,3);

    // enthalpy_extrapolated_jac+=cell.dExtrapolateQuant(Physquant::enthalpyFlow);
    enthalpy_extrapolated_jac+=JacCol(cell.mHbc(),-eye);
    enthalpy_extrapolated_jac+=JacCol(cell.mbc(),fDFT*diagmat(cp(cell)*cell.Tbc().tdata())*iDFT);
    enthalpy_extrapolated_jac+=JacCol(cell.Tbc(),fDFT*diagmat(cp(cell)*cell.mbc().tdata())*iDFT);
    jac+=enthalpy_extrapolated_jac;

  }
  void PressureBc::show(us detailnr) const {
    TRACE(5,"PressureBc::show()");
    checkInit();
    const char* side;
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



