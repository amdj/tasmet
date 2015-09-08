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

  // The functions adiabaticTemp() and coldTemp are defined in varutils.h

  // If only pressure is given, solid temperature is set to gc->T0 and
  // fluid temperature to adiabatic compression/expansion
  PressureBc::PressureBc(const string& segid,Pos position,const var& pres):
    PressureBc(segid,position,pres,coldTemp(pres))
  {}				// 
  // If no fluid temperature is given, the fluid temperature follows
  // adiabatic compression/expansion
  PressureBc::PressureBc(const string& segid,Pos position,const var& pres,\
			 const var& stemp):
    PressureBc(segid,position,pres,stemp,adiabaticTemp(pres))
  {}
  // The full constructor
  PressureBc::PressureBc(const string& segid,Pos position,const var& pres,\
			 const var& stemp,const var& temp):
    TubeBc(segid,position),
    p_prescribed(pres),
    prescribeT(temp),
    prescribeTs(stemp)    
  {
    TRACE(8,"PressureBc full constructor");
  }

  PressureBc::PressureBc(const PressureBc& other,const TaSystem& sys):
    TubeBc(other,sys),
    p_prescribed(other.p_prescribed),
    prescribeT(other.prescribeT),
    prescribeTs(other.prescribeTs)
  {
    TRACE(8,"PressureBc copy constructor");
    assert(gc);

    // Decouple from old globalconf pointer
    p_prescribed.setGc(*gc);
    prescribeT.setGc(*gc);
    prescribeTs.setGc(*gc);
  }
 
  void PressureBc::updateNf(){
    p_prescribed.updateNf();
    prescribeT.updateNf();
    prescribeTs.updateNf();
  }
  void PressureBc::setEqNrs(us firsteqnr){
    TRACE(2,"Pressure::setEqNrs()");
    TubeBc::setEqNrs(firsteqnr);
    const BcCell& cell=t->bcCell(pos);
    prescribeT.set(firsteqnr+2*Ns,cell.Tbc());
    prescribeTs.set(firsteqnr+3*Ns,cell.Tsbc());
  }
  us PressureBc::getNEqs() const {
    TRACE(15,"PressureBc::getNEqs()");
    if(t->hasSolid())
      return 4*Ns;
    else
      return 3*Ns;
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
      errorM+=Wddt*(t->dragResistance().drag(cell));
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
      errorM+=Wddt*(t->dragResistance().drag(cell));
      #endif
    }
    // Momentum equation
    error.subvec(0,Ns-1)=errorM;

    // Extrapolation of enthalpy flow
    error.subvec(Ns,2*Ns-1)=-cell.mHbc()();
    error.subvec(Ns,2*Ns-1)+=cell.extrapolateQuant(Varnr::mH);

    // Prescribed temperature
    error.subvec(2*Ns,3*Ns-1)=prescribeT.error();
    if(t->hasSolid())
      error.subvec(3*Ns,4*Ns-1)=prescribeTs.error();    
    return error;
  }
  void PressureBc::jac(Jacobian& jac) const{
    TRACE(8,"PressureBc::jac()");

    const BcCell& cell=t->bcCell(pos);
    // Momentum equation Jacobian
    if(pos==Pos::left)    {
      d Wddt=cell.vx;
      JacRow jacr(firsteqnr,4);
      jacr+=JacCol(cell.mbc(),Wddt*DDTfd);
      jacr+=JacCol(cell.p(),cell.vSf*eye);
      jacr+=JacCol(cell.mu(),eye);
      jacr+=-cell.dExtrapolateQuant(Varnr::mu);
      #ifndef NODRAG
      jacr+=JacCol(cell.ml(),Wddt*(t->dragResistance().dm(cell)));
      #endif  // NODRAG
      jac+=jacr;
    }
    else{
      d Wddt=cell.xr-cell.vx;
      JacRow jacr(firsteqnr,4);
      jacr+=JacCol(cell.mbc(),Wddt*DDTfd);
      jacr+=JacCol(cell.p(),-cell.vSf*eye);
      jacr+=JacCol(cell.mu(),-eye);
      jacr+=cell.dExtrapolateQuant(Varnr::mu);
      #ifndef NODRAG
      jacr+=JacCol(cell.mr(),Wddt*(t->dragResistance().dm(cell)));
      #endif  // NODRAG
      jac+=jacr;
    }
    // Prescribed enthalpy flow.
    JacRow enthalpy_extrapolated_jac(firsteqnr+Ns,3);
    enthalpy_extrapolated_jac+=cell.dExtrapolateQuant(Varnr::mH);
    enthalpy_extrapolated_jac+=JacCol(cell.mHbc(),-eye);
    // enthalpy_extrapolated_jac+=JacCol(cell.mbc(),fDFT*diagmat(cp(cell)*cell.Tbc().tdata())*iDFT);
    // enthalpy_extrapolated_jac+=JacCol(cell.Tbc(),fDFT*diagmat(cp(cell)*cell.mbc().tdata())*iDFT);
    jac+=enthalpy_extrapolated_jac;

    // Prescribed temperature Jacobian part
    jac+=prescribeT.jac();
    // Prescribed solid temperature Jacobian part
    if(t->hasSolid())
      jac+=prescribeTs.jac();

  }
  void PressureBc::show(us detailnr) const {
    TRACE(5,"PressureBc::show()");
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



