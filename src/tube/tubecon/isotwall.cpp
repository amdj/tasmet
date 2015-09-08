#include "isotwall.h"
#include "jacobian.h"
#include "bccell.h"
#include "tube.h"

#define iDFT (gc->iDFT)
#define fDFT (gc->fDFT)
#define DDTfd (gc->DDTfd)
#define Ns (gc->Ns())


namespace tube{

  using tasystem::var;
  using tasystem::TaSystem;
  using tasystem::Globalconf;
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  IsoTWall::IsoTWall(const string& segid,Pos position,const var& Tbc,bool arbitrateMass):
    TubeBc(segid,position),
    Tbc(Tbc),
    Tsbc(Tbc),
    arbitrateMass(arbitrateMass)
  {
    TRACE(15,"IsoTWall::IsoTWall()");
  }
  IsoTWall::IsoTWall(const IsoTWall& o,const TaSystem& sys):
    TubeBc(o,sys),
    massflowzero(o.massflowzero),
    enthalpyflowzero(o.enthalpyflowzero),
    Tbc(o.Tbc),
    Tsbc(o.Tsbc),
    arbitrateMass(o.arbitrateMass)
  {
    TRACE(15,"IsoTWall::IsoTWall(copy)");
    massflowzero.setGc(*gc);
    enthalpyflowzero.setGc(*gc);
    Tbc.setGc(*gc);
    Tsbc.setGc(*gc);
  }
  void IsoTWall::updateNf(){
    TRACE(15,"IsoTWall::updateNf()");
    massflowzero.updateNf();
    enthalpyflowzero.updateNf();
    Tbc.updateNf();
    Tsbc.updateNf();
  }
  int IsoTWall::arbitrateMassEq() const {
    TRACE(15,"AdiabaticWall::arbitrateMassEq()");
    if(arbitrateMass)
      return firsteqnr;
    else
      return -1;
  }
  void IsoTWall::show(us i) const {
    TRACE(5,"IsoTWall::show()");
    cout << "IsoTWall boundary condition set at "<<posWord(pos) <<" side of segment .\n";
  }
  us IsoTWall::getNEqs() const {
    if(t->hasSolid())
      return 4*Ns;
    else
      return 3*Ns;
  }
  void IsoTWall::setEqNrs(us firsteqnr){
    TRACE(2,"IsoTWall::setEqNrs()");
    TubeBc::setEqNrs(firsteqnr);
    var zero(*gc,0);             // zero at all the time
    const BcCell& cell=t->bcCell(pos);
    massflowzero.set(firsteqnr,cell.mbc(),zero);
    enthalpyflowzero.set(firsteqnr+Ns,cell.mHbc(),zero);
    Tbc.set(firsteqnr+2*Ns,cell.Tbc()); 
    if(t->hasSolid())
      Tsbc.set(firsteqnr+3*Ns,cell.Tsbc()); 
  }
  vd IsoTWall::error() const {
    TRACE(15,"IsoTWall::error()");
    vd error(getNEqs());

    // Add all individual error parts to vector
    error.subvec(0,Ns-1)=massflowzero.error();
    error.subvec(Ns,2*Ns-1)=enthalpyflowzero.error();    
    error.subvec(2*Ns,3*Ns-1)=Tbc.error();
    if(t->hasSolid())
      error.subvec(3*Ns,4*Ns-1)=Tsbc.error();
    return error;
  }
  void IsoTWall::jac(Jacobian& jac) const {
    TRACE(15,"IsoTWall::jac()");

    jac+=massflowzero.jac();
    jac+=enthalpyflowzero.jac();
    jac+=Tbc.jac();
    if(t->hasSolid())
      jac+=Tsbc.jac();
  }
} // namespace tube

