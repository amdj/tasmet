// ductpistonconnector.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "ductpistonconnector.h"
#include "tasystem.h"
#include "piston.h"
#include "staticmsg.h"
#include "duct.h"
#include "bccell.h"
#include "jacobian.h"
#include "jacrow.h"

#define fDFT (gc->fDFT())
#define iDFT (gc->iDFT())
#define DDTfd (gc->DDTfd())

#define Ns (gc->Ns())
#define eye (arma::eye(Ns,Ns))

namespace{
    common::StaticMsg<> msg;
}
namespace mech {
  using tasystem::JacRow;
  using tasystem::JacCol;
  using tasystem::Jacobian;
  using tasystem::var;
  using tasystem::Globalconf;


  
  DuctPistonConnector::DuctPistonConnector(const string& ductid,Pos ductpos,\
                                           const string& pistonid,Pos pistonpos,\
                                           d KDuctPiston,d KPistonDuct):
    Connector(),
    pistonid(pistonid),
    ductid(ductid),
    pistonPos(pistonpos),
    ductPos(ductpos),
    KDuctPiston(KDuctPiston),
    KPistonDuct(KPistonDuct)
  {
    TRACE(15,"DuctPistonConnector::DuctPistonConnector()");
    if(!(KDuctPiston>=0) || !(KDuctPiston<=1))
       throw MyError(msg("Illegal value given for KDuctPiston. Should be >=0"
                         " and <=1. Value given: %0.2e",KDuctPiston));
    if(!(KPistonDuct>=0) || !(KPistonDuct<=1))
       throw MyError(msg("Illegal value given for KPistonDuct. Should be >=0"
                         " and <=1. Value given: %0.2e",KPistonDuct));
  }
  DuctPistonConnector::DuctPistonConnector(const DuctPistonConnector& other,const tasystem::TaSystem& sys):
    Connector(other,sys),
    pistonid(other.pistonid),
    ductid(other.ductid),
    pistonPos(other.pistonPos),
    ductPos(other.ductPos),
    KDuctPiston(other.KDuctPiston),
    KPistonDuct(other.KPistonDuct)
  {
    TRACE(15,"DuctPistonConnector::DuctPistonConnector()");
    if((piston=&sys.getPiston(pistonid))==nullptr){
      throw MyError(msg("Segment %s is not of type Piston",pistonid.c_str()));
    }
    if((duct=&sys.getDuct(ductid))==nullptr){
      throw MyError(msg("Segment %g is not of type Duct",ductid.c_str()));
    }
    
    // Check if piston is already connected. If so, throw error
    if(piston->isConnected(pistonPos))
       throw MyError(msg("Error: piston %d is already connected to "
                         "a segment at side %s.",pistonid.c_str(),duct::posWord(pistonPos)));

    // Tell the Piston that it is connected at some side
    piston->setConnected(pistonPos);
  }

  vd DuctPistonConnector::error() const {
    TRACE(15,"vd DuctPistonConnector::error()");
    const duct::BcCell& ductcell=duct->bcCell(ductPos);
    const PistonConfiguration& pc=piston->getPc();
    // Important: the sum of flows OUT of each segment should be
    // zero. However, if the position is left, we count out flow as
    // the negative of inflow
    int signduct=ductPos==Pos::left?-1:1;


    vd error(getNEqs());
    // ***** First equation: mass balance
    error.subvec(0,Ns-1)=piston->m(pistonPos)()\
      +signduct*ductcell.mbc()();

    // ***** Second equation, extrapolated mH equals 
    // ***** mH
    error.subvec(Ns,2*Ns-1)=ductcell.extrapolateQuant(Varnr::mH)
      -ductcell.mHbc()();

    // ***** Third equation: energy balance
     error.subvec(2*Ns,3*Ns-1)=signduct*ductcell.mHbc()()+\
      signduct*ductcell.extrapolateQuant(Varnr::Q)\
      +piston->mH(pistonPos)();    

    // ***** Fourth equation: conduction to (from) piston equals
    // ***** conduction from duct

   // Compute conduction path length
    d Lcon=pistonPos==Pos::left?pc.V0l/pc.Sl:pc.V0r/pc.Sr;
    // Compute cross sectional area of duct exit
    d Sf=ductcell.Sfbc();
    // Temperature of gas in piston volume in time domain
    const vd& pTt=piston->T(pistonPos).tdata();
    // Temperature of gas at duct boundary in time domain
    const vd& tTt=ductcell.Tbc().tdata();
    // Kappa of piston in time domain.
    vd kappat=gc->gas().kappa(pTt);


    vd Qpistontoduct=fDFT*((Sf/Lcon)*kappat%\
                           (pTt-tTt));

    error.subvec(3*Ns,4*Ns-1)=Qpistontoduct\
      -signduct*ductcell.extrapolateQuant(Varnr::Q);


    // ***** Fifth equation: something with pressure
    error.subvec(4*Ns,5*Ns-1)=ductcell.pbc()()	\
      -piston->p(pistonPos)();
    return error;
  }
  void DuctPistonConnector::jac(tasystem::Jacobian& jac) const{
    TRACE(15,"void DuctPistonConnector::jac()");
    assert(piston&&duct);
    assert(piston&&duct);
    const duct::BcCell& ductcell=duct->bcCell(ductPos);
    const PistonConfiguration& pc=piston->getPc();
    // Important: the sum of flows OUT of each segment should be
    // zero. However, if the position is left, we count out flow as
    // the negative of inflow
    int signduct=ductPos==Pos::left?-1:1;

    // Temperature of gas in piston volume in time domain
    const vd& pTt=piston->T(pistonPos).tdata();
    // Temperature of gas at duct boundary in time domain
    const vd& tTt=ductcell.Tbc().tdata();
    // Mass flow out of piston in time domain
    const vd& pmt=piston->m(pistonPos).tdata();
    // Kappa of piston in time domain.
    vd kappat=gc->gas().kappa(pTt);


    // First equation: mass balance
    JacRow massjac(firsteqnr,2);
    massjac+=JacCol(piston->m(pistonPos),eye);
    massjac+=JacCol(ductcell.mbc(),signduct*eye);
    jac+=massjac;

    // ***** Second equation, extrapolated mH equals 
    // ***** mH
    JacRow enthexpjac(firsteqnr+Ns,4);
    enthexpjac+=ductcell.dExtrapolateQuant(Varnr::mH);
    enthexpjac+=JacCol(ductcell.mHbc(),-eye);
    jac+=enthexpjac;

    // ***** Third equation: energy balance
    JacRow enbalancejac(firsteqnr+2*Ns,6);
    enbalancejac+=JacCol(ductcell.mHbc(),signduct*eye);
    enbalancejac+=(ductcell.dExtrapolateQuant(Varnr::Q)*=signduct);
    enbalancejac+=JacCol(piston->mH(pistonPos),eye);
    jac+=enbalancejac;

    // ***** Fourth equation: conduction to (from) piston equals
    // ***** conduction from duct

    // Compute conduction path length
    d Lcon=(pistonPos==Pos::left)?pc.V0l/pc.Sl:pc.V0r/pc.Sr;
    // Compute cross sectional area of duct exit
    d Sf=ductcell.Sfbc();

    JacRow QisQjac(firsteqnr+3*Ns,3);
    QisQjac+=JacCol(piston->T(pistonPos),                       \
                         fDFT*(Sf/Lcon)*diagmat(kappat)*iDFT);
    QisQjac+=JacCol(ductcell.Tbc(),                             \
                         -fDFT*(Sf/Lcon)*diagmat(kappat)*iDFT);
    QisQjac+=(ductcell.dExtrapolateQuant(Varnr::Q)*=-signduct);
    jac+=QisQjac;

    // ***** Fifth equation: something with pressure
    JacRow pispjac(firsteqnr+4*Ns,3);
    pispjac+=JacCol(ductcell.pbc(),eye);
    pispjac+=JacCol(piston->p(pistonPos),-eye);
    jac+=pispjac;
  }
  void DuctPistonConnector::setEqNrs(us firsteqnr){
    TRACE(15,"vd DuctPistonConnector::setEqNrs()");
    this->firsteqnr=firsteqnr;
  }
  us DuctPistonConnector::getNEqs() const{
    TRACE(15,"us DuctPistonConnector::getNEqs()");
    return 5*Ns;
  }
  void DuctPistonConnector::updateNf() {
    TRACE(15,"void DuctPistonConnector::updateNf()");
  }
  void DuctPistonConnector::show(us detailnr) const{

    std::cout << "DuctPistonConnector which connects duct " << ductid <<
      " at the " << duct::posWord(ductPos) << " side to Piston " << pistonid <<
      " on the " << duct::posWord(pistonPos) << " side.\n";
  }

} // namespace mech
//////////////////////////////////////////////////////////////////////

