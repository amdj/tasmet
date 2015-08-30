#include "pressurebc.h"
#include "tasystem.h"
#include "jacobian.h"
#include "bccell.h"
#include "constants.h"
#include "state.h"
#include "tube.h"

#define fDFT (gc->fDFT)
#define iDFT (gc->iDFT)
#define DDTfd (gc->DDTfd)

#define Ns (gc->Ns())
#define eye (eye(Ns,Ns))
#include "varutils.h"

namespace tube{
  using tasystem::var;
  using tasystem::coldTemp;
  using tasystem::adiabaticTemp;
  using tasystem::TaSystem;
  using tasystem::Globalconf;
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  PressureBc::PressureBc(const string& segid,Pos position,const var& pres,const var& temp,const var& stemp):
    TubeBc(segid,position),
    p_prescribed(pres),
    prescribeT(temp)
  {
    
    TRACE(8,"PressureBc full constructor");

  }
  PressureBc::PressureBc(const string& segid,Pos position,const var& pres,const var& temp):
    PressureBc(segid,position,pres,temp,coldTemp(pres))
  {}
  PressureBc::PressureBc(const string& segid,Pos position,const var& pres):
    PressureBc(segid,position,pres,adiabaticTemp(pres))
  {}
  PressureBc::PressureBc(const PressureBc& other,const TaSystem& sys):
    TubeBc(other,sys),
    p_prescribed(other.p_prescribed),
    prescribeT(other.prescribeT)
  {
    TRACE(8,"PressureBc copy constructor");
    assert(gc);

    // Decouple from old globalconf pointer
    p_prescribed.setGc(*gc);
    prescribeT.setGc(*gc);
    setInit(true);
  }
 
  void PressureBc::updateNf(){
    p_prescribed.updateNf();
    prescribeT.updateNf();
  }
  void PressureBc::setEqNrs(us firsteqnr){
    TRACE(2,"Pressure::setEqNrs()");
    TubeBc::setEqNrs(firsteqnr);
    const BcCell& cell=t->bcCell(pos);
    prescribeT.set(firsteqnr+Ns,cell.Tbc());
  }

  vd PressureBc::error() const {
    TRACE(15,"PressureBc::error()");
    vd error(getNEqs());
    const BcCell& cell=t->bcCell(pos);

    vd errorM(Ns,fillwith::zeros);
    if(pos==Pos::left)    {
      d Wddt=cell.vx;
      errorM+=Wddt*DDTfd*cell.mbc()();
      errorM+=cell.vSf*cell.p()();
      errorM-=cell.Sfl*p_prescribed();
      errorM+=cell.mu()();       // 
      errorM-=cell.extrapolateQuant(Varnr::mu);
      #ifndef NODRAG
      errorM+=Wddt*(t->getDragResistance().drag(cell));
      #endif
    }
    else{
      d Wddt=cell.xr-cell.vx;
      errorM+=Wddt*DDTfd*cell.mbc()();
      errorM-=cell.vSf*cell.p()();
      errorM+=cell.Sfr*p_prescribed();
      errorM-=cell.mu()();
      errorM+=cell.extrapolateQuant(Varnr::mu);      
      #ifndef NODRAG
      errorM+=Wddt*(t->getDragResistance().drag(cell));
      #endif
    }      
    error.subvec(0,Ns-1)=errorM;
    error.subvec(Ns,2*Ns-1)=prescribeT.error();

    error.subvec(2*Ns,3*Ns-1)=-cell.mHbc()();
    error.subvec(2*Ns,3*Ns-1)+=cell.extrapolateQuant(Varnr::mH);
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
      JacRow jacr(firsteqnr,4);
      jacr+=JacCol(cell.mbc(),Wddt*DDTfd);
      jacr+=JacCol(cell.p(),cell.vSf*eye);
      jacr+=JacCol(cell.mu(),eye);
      jacr+=-cell.dExtrapolateQuant(Varnr::mu);
      #ifndef NODRAG
      jacr+=JacCol(cell.ml(),Wddt*(t->getDragResistance().dm(cell)));
      #endif  // NODRAG
      jac+=jacr;
    }
    else{
      d Wddt=cell.xr-cell.vx;
      // VARTRACE(25,Wddt);
      JacRow jacr(firsteqnr,4);
      jacr+=JacCol(cell.mbc(),Wddt*DDTfd);
      jacr+=JacCol(cell.p(),-cell.vSf*eye);
      jacr+=JacCol(cell.mu(),-eye);
      jacr+=cell.dExtrapolateQuant(Varnr::mu);
      #ifndef NODRAG
      jacr+=JacCol(cell.mr(),Wddt*(t->getDragResistance().dm(cell)));
      #endif  // NODRAG
      jac+=jacr;
    }
    // Prescribed enthalpy flow.
    JacRow enthalpy_extrapolated_jac(firsteqnr+2*Ns,3);

    enthalpy_extrapolated_jac+=cell.dExtrapolateQuant(Varnr::mH);
    enthalpy_extrapolated_jac+=JacCol(cell.mHbc(),-eye);
    // enthalpy_extrapolated_jac+=JacCol(cell.mbc(),fDFT*diagmat(cp(cell)*cell.Tbc().tdata())*iDFT);
    // enthalpy_extrapolated_jac+=JacCol(cell.Tbc(),fDFT*diagmat(cp(cell)*cell.mbc().tdata())*iDFT);
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
    cout << "PressureBc boundary condition set at "<<side <<" side of segment "<< segid<<".\n";
    if(detailnr>1){
      cout << "Prescribed pressure:" << p_prescribed() << "\n";
      cout << "Prescribed fluid temperature:" << prescribeT.getVals()() << "\n";      
      // cout << "Prescribed solid temperature:" << prescribeTs.getVals()() << "\n";      
    }
  } // show

} // namespace tube



