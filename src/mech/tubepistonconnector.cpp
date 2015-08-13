// tubepistonconnector.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "tubepistonconnector.h"
#include "tasystem.h"
#include "piston.h"
#include "staticmsg.h"
#include "tube.h"
#include "bccell.h"
#include "jacobian.h"
#include "jacrow.h"

#define fDFT (gc->fDFT)
#define iDFT (gc->iDFT)
#define DDTfd (gc->DDTfd)

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


  
  TubePistonConnector::TubePistonConnector(us tubenr,Pos tubepos,\
                                           us pistonnr,Pos pistonpos,\
                                           d KTubePiston,d KPistonTube):
    Connector(),
    pistonNr(pistonnr),
    tubeNr(tubenr),
    pistonPos(pistonpos),
    tubePos(tubepos),
    KTubePiston(KTubePiston),
    KPistonTube(KPistonTube)
  {
    TRACE(15,"TubePistonConnector::TubePistonConnector()");
    if(!(KTubePiston>=0) || !(KTubePiston<=1))
       throw MyError(msg("Illegal value given for KTubePiston. Should be >=0"
                         " and <=1. Value given: %0.2e",KTubePiston));
    if(!(KPistonTube>=0) || !(KPistonTube<=1))
       throw MyError(msg("Illegal value given for KPistonTube. Should be >=0"
                         " and <=1. Value given: %0.2e",KPistonTube));
  }
  TubePistonConnector::TubePistonConnector(const TubePistonConnector& other,const tasystem::TaSystem& sys):
    Connector(other,sys),
    pistonNr(other.pistonNr),
    tubeNr(other.tubeNr),
    pistonPos(other.pistonPos),
    tubePos(other.tubePos),
    KTubePiston(other.KTubePiston),
    KPistonTube(other.KPistonTube)
  {
    TRACE(15,"TubePistonConnector::TubePistonConnector()");
    if((piston=&sys.getPiston(pistonNr))==nullptr){
      throw MyError(msg("Segment %g is not of type Piston",pistonNr));
    }
    if((tube=&sys.getTube(tubeNr))==nullptr){
      throw MyError(msg("Segment %g is not of type Tube",tubeNr));
    }
    
    // Check if piston is already connected. If so, throw error
    if(piston->isConnected(pistonPos))
       throw MyError(msg("Error: piston %d is already connected to "
                         "a segment at side %s.",pistonNr,tube::posWord(pistonPos)));

    // Tell the Piston that it is connected at some side
    piston->setConnected(pistonPos);
  }

  vd TubePistonConnector::error() const {
    TRACE(15,"vd TubePistonConnector::error()");
    const tube::BcCell& tubecell=tube->bcCell(tubePos);
    const PistonConfiguration& pc=piston->getPc();
    // Important: the sum of flows OUT of each segment should be
    // zero. However, if the position is left, we count out flow as
    // the negative of inflow
    int signtube=tubePos==Pos::left?-1:1;


    vd error(getNEqs());
    // ***** First equation: mass balance
    error.subvec(0,Ns-1)=piston->m(pistonPos)()\
      +signtube*tubecell.mbc()();

    // ***** Second equation, extrapolated mH equals 
    // ***** mH
    error.subvec(Ns,2*Ns-1)=tubecell.extrapolateQuant(Varnr::mH)
      -tubecell.mHbc()();

    // ***** Third equation: energy balance
     error.subvec(2*Ns,3*Ns-1)=signtube*tubecell.mHbc()()+\
      signtube*tubecell.extrapolateQuant(Varnr::Q)\
      +piston->mH(pistonPos)();    

    // ***** Fourth equation: conduction to (from) piston equals
    // ***** conduction from tube

   // Compute conduction path length
    d Lcon=pistonPos==Pos::left?pc.V0l/pc.Sl:pc.V0r/pc.Sr;
    // Compute cross sectional area of tube exit
    d Sf=tubecell.Sfbc();
    // Temperature of gas in piston volume in time domain
    const vd& pTt=piston->T(pistonPos).tdata();
    // Temperature of gas at tube boundary in time domain
    const vd& tTt=tubecell.Tbc().tdata();
    // Kappa of piston in time domain.
    vd kappat=gc->gas().kappa(pTt);


    vd Qpistontotube=fDFT*((Sf/Lcon)*kappat%\
                           (pTt-tTt));

    error.subvec(3*Ns,4*Ns-1)=Qpistontotube\
      -signtube*tubecell.extrapolateQuant(Varnr::Q);


    // ***** Fifth equation: something with pressure
    error.subvec(4*Ns,5*Ns-1)=tubecell.extrapolateQuant(Varnr::p)\
      -piston->p(pistonPos)();
    return error;
  }
  void TubePistonConnector::jac(tasystem::Jacobian& jac) const{
    TRACE(15,"void TubePistonConnector::jac()");
    assert(piston&&tube);
    assert(piston&&tube);
    const tube::BcCell& tubecell=tube->bcCell(tubePos);
    const PistonConfiguration& pc=piston->getPc();
    // Important: the sum of flows OUT of each segment should be
    // zero. However, if the position is left, we count out flow as
    // the negative of inflow
    int signtube=tubePos==Pos::left?-1:1;

    // Temperature of gas in piston volume in time domain
    const vd& pTt=piston->T(pistonPos).tdata();
    // Temperature of gas at tube boundary in time domain
    const vd& tTt=tubecell.Tbc().tdata();
    // Mass flow out of piston in time domain
    const vd& pmt=piston->m(pistonPos).tdata();
    // Kappa of piston in time domain.
    vd kappat=gc->gas().kappa(pTt);


    // First equation: mass balance
    JacRow massjac(firsteqnr,2);
    massjac+=JacCol(piston->m(pistonPos),eye);
    massjac+=JacCol(tubecell.mbc(),signtube*eye);
    jac+=massjac;

    // ***** Second equation, extrapolated mH equals 
    // ***** mH
    JacRow enthexpjac(firsteqnr+Ns,4);
    enthexpjac+=tubecell.dExtrapolateQuant(Varnr::mH);
    enthexpjac+=JacCol(tubecell.mHbc(),-eye);
    jac+=enthexpjac;

    // ***** Third equation: energy balance
    JacRow enbalancejac(firsteqnr+2*Ns,6);
    enbalancejac+=JacCol(tubecell.mHbc(),signtube*eye);
    enbalancejac+=(tubecell.dExtrapolateQuant(Varnr::Q)*=signtube);
    enbalancejac+=JacCol(piston->mH(pistonPos),eye);
    jac+=enbalancejac;

    // ***** Fourth equation: conduction to (from) piston equals
    // ***** conduction from tube

    // Compute conduction path length
    d Lcon=(pistonPos==Pos::left)?pc.V0l/pc.Sl:pc.V0r/pc.Sr;
    // Compute cross sectional area of tube exit
    d Sf=tubecell.Sfbc();

    JacRow QisQjac(firsteqnr+3*Ns,3);
    QisQjac+=JacCol(piston->T(pistonPos),                       \
                         fDFT*(Sf/Lcon)*diagmat(kappat)*iDFT);
    QisQjac+=JacCol(tubecell.Tbc(),                             \
                         -fDFT*(Sf/Lcon)*diagmat(kappat)*iDFT);
    QisQjac+=(tubecell.dExtrapolateQuant(Varnr::Q)*=-signtube);
    jac+=QisQjac;

    // ***** Fifth equation: something with pressure
    JacRow pispjac(firsteqnr+4*Ns,3);
    pispjac+=tubecell.dExtrapolateQuant(Varnr::p);
    pispjac+=JacCol(piston->p(pistonPos),-eye);
    jac+=pispjac;
  }
  void TubePistonConnector::setEqNrs(us firsteqnr){
    TRACE(15,"vd TubePistonConnector::setEqNrs()");
    this->firsteqnr=firsteqnr;
  }
  us TubePistonConnector::getNEqs() const{
    TRACE(35,"us TubePistonConnector::getNEqs()");
    return 5*Ns;
  }
  void TubePistonConnector::updateNf() {
    TRACE(15,"void TubePistonConnector::updateNf()");
  }
  void TubePistonConnector::show(us detailnr) const{

    std::cout << "TubePistonConnector which connects tube " << tubeNr <<
      " at the " << tube::posWord(tubePos) << " side to Piston " << pistonNr <<
      " on the " << tube::posWord(pistonPos) << " side.\n";
  }

} // namespace mech
//////////////////////////////////////////////////////////////////////

