// tubeconnector.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
// Connect two tube-type segments by using conservation of mass,
// momentum and energy. Number of equations 
//////////////////////////////////////////////////////////////////////
// #define TRACERPLUS (20)
#include "tubeconnector.h"
#include "tasystem.h"
#include "tube.h"
#include "bccell.h"
#include "constants.h"
#include "jacobian.h"
#include "energy.h"

#define fDFT (gc->fDFT)
#define iDFT (gc->iDFT)
#define DDTfd (gc->DDTfd)

#define Ns (gc->Ns())
#define eye (arma::eye(Ns,Ns))

namespace tube {
  using tasystem::TaSystem;
  using tasystem::JacRow;
  using tasystem::JacCol;
  // Number of equations corresponding to this connection:
  const int Neq=6;

  SimpleTubeConnector::SimpleTubeConnector(us seg1,Pos pos1,\
                                           us seg2,Pos pos2,d K1to2,d K2to1):
    K1to2(K1to2),
    K2to1(K2to1)
  {
    TRACE(15,"SimpleTubeConnector::SimpleTubeConnector()");
  
    if(max(seg1,seg2)>constants::maxsegs)
      throw MyError("Too high segment number given");
    if(seg1==seg2)
      throw MyError("Segments cannot be the same");
    segnrs[0]=seg1;
    segnrs[1]=seg2;
    pos[0]=pos1;
    pos[1]=pos2;
      
  }
  SimpleTubeConnector::SimpleTubeConnector(const SimpleTubeConnector& o,
                                           const TaSystem& sys):
    Connector(o,sys),
    segnrs(o.segnrs),
    pos(o.pos),
    K1to2(o.K1to2),
    K2to1(o.K2to1)
  {
    
    bccells[0]=&sys.getTube(segnrs[0]).bcCell(pos[0]);
    bccells[1]=&sys.getTube(segnrs[1]).bcCell(pos[1]);
    assert(bccells[0]&&bccells[1]);

    if(pos[0]==Pos::left) 
      out[0]=-1;
    if(pos[1]==Pos::left)
      out[1]=-1;
    
    setInit(true);
  }

  #define Defs						\
    d T0=gc->T0();					\
    d gamma=gc->gas().gamma(T0);			\
    d p0=gc->p0();					\
    d cp=gc->gas().cp(T0);				\
    d Rs=gc->gas().Rs();				\
    d Sf0=bccells[0]->Sfbc();				\
    d Sf1=bccells[1]->Sfbc();				\
    const vd p0tbc=bccells[0]->pbc().tdata()+p0;	\
    const vd p1tbc=bccells[1]->pbc().tdata()+p0;	\
    const vd& T0tbc=bccells[0]->Tbc().tdata();		\
    const vd& T1tbc=bccells[1]->Tbc().tdata();		\
    const vd& m0tbc=bccells[0]->mbc().tdata();		\
    const vd& m1tbc=bccells[1]->mbc().tdata();		\
    const vd& T0t=bccells[0]->T().tdata();		\
    const vd& T1t=bccells[1]->T().tdata();		\
    const vd& rho0tbc=bccells[0]->rhobc().tdata();	\
    const vd& rho1tbc=bccells[1]->rhobc().tdata();	\
  const vd& u0tbc=bccells[0]->ubc().tdata();		\
  const vd& u1tbc=bccells[1]->ubc().tdata();		

  vd SimpleTubeConnector::error() const{
    TRACE(10,"SimpleTubeConnector::error()");

    us nr=0;
    vd error(Neq*Ns);

    Defs
    
      // First equation: mass balance
      TRACE(15," First equation: mass balance");
    error.subvec(Ns*nr,Ns*(nr+1)-1)=\
      out[0]*bccells[0]->mbc()()\
      +out[1]*bccells[1]->mbc()();
    nr++;

    // Second equation: enthalpy flow balance
    TRACE(15," Second equation: enthalpy flow balance");
    error.subvec(Ns*nr,Ns*(nr+1)-1)=\
      out[0]*bccells[0]->mHbc()()\
      +out[1]*bccells[1]->mHbc()();
    nr++;

    // Third equation: Enthalpy flow at the interface equals the average of the
    // enthalpy flows of both tubes
    TRACE(15," Third equation: Enthalpy flow at the interface equals the average of the");
    error.subvec(Ns*nr,Ns*(nr+1)-1)=\
      0.5*out[0]*bccells[0]->extrapolateQuant(Varnr::mH)\
      -0.5*out[1]*bccells[1]->extrapolateQuant(Varnr::mH)\
      -out[0]*bccells[0]->mHbc()();
    nr++;

    // Fourth equation: continuity of total enthalpy
    error.subvec(Ns*nr,Ns*(nr+1)-1)=\
      fDFT*(cp*T0tbc+0.5*pow(u0tbc,2)-cp*T1tbc-0.5*pow(u1tbc,2));
      // fDFT*(kappaSft%(T0t-T1t)/dx)\
      // -out[0]*bccells[0]->extrapolateQuant(Varnr::Q);
    nr++;      

    // Fifth equation: heat out of one segment goes into the other
    TRACE(15," Fifth equation: heat out of one segment goes into the other");
    error.subvec(Ns*nr,Ns*(nr+1)-1)=\
      out[0]*bccells[0]->extrapolateQuant(Varnr::Q)\
      +out[1]*bccells[1]->extrapolateQuant(Varnr::Q);
    nr++;

    // 6th equation: Change in exergy due to minor loss
    TRACE(15," 6th equation: Change in exergy due to minor loss");
    d one_eight=1.0/8.0;
    // vd DeltaEx=fDFT*log((p1tbc/p0tbc)%pow(rho0tbc/rho1tbc,gamma));
 
    vd DeltaEx=fDFT*log((p0tbc/p1tbc)%pow(T1tbc/T0tbc,gamma/(gamma-1)));
    // vd DeltaEx=fDFT*log(p1tbc/p0tbc);
    // vd DeltaEx=fDFT*(p1tbc-p0tbc);
    VARTRACE(15,DeltaEx);
    // vd u0t=(iDFT*bccells[0]->extrapolateQuant(Varnr::mu))/bccells[0]->mbc().tdata();
    // vd u1t=(iDFT*bccells[1]->extrapolateQuant(Varnr::mu))/bccells[1]->mbc().tdata();
    // vd minus_minorLoss=K1to2*(T0/T1t)*one_eight	\
    //   *pow(out[0]*abs(u1t)+u1t,2)-\
    //   K2to1*(T0/T2t)*one_eight\
    //   *pow(out[1]*abs(u2t)+u2t,2);
    // vd minus_minorLoss=zeros(Ns);
    error.subvec(Ns*nr,Ns*(nr+1)-1)=\
      DeltaEx;//		    \
    //+minus_minorLoss;
    // error.subvec(Ns*nr,Ns*(nr+1)-1)=\
    // Sfgem*(bccells[1]->extrapolateQuant(Varnr::p)\
    // 	     -bccells[0]->extrapolateQuant(Varnr::p));//+
    // bccells[1]->extrapolateQuant(Varnr::mu)
    // -bccells[0]->extrapolateQuant(Varnr::mu);
    nr++;

    return error;
  }
  void SimpleTubeConnector::jac(tasystem::Jacobian& jac) const {
    TRACE(15,"SimpleTubeconnector::jac()");
    us eqnr=firsteqnr;

    Defs
    
      // First equation: mass balance
      JacRow mjac(eqnr,2);
    mjac+=JacCol(bccells[0]->mbc(),out[0]*eye);
    mjac+=JacCol(bccells[1]->mbc(),out[1]*eye);
    jac+=mjac;
    eqnr+=Ns;

    // Second equation: enthalpy flow balance
    JacRow mHjac(eqnr,2);
    mHjac+=JacCol(bccells[0]->mHbc(),out[0]*eye);
    mHjac+=JacCol(bccells[1]->mHbc(),out[1]*eye);
    jac+=mHjac;
    eqnr+=Ns;

    // // Third equation: Enthalpy flow at the interface equals the average of the
    // // enthalpy flows of both tubes
    JacRow mHjacint(eqnr,5);
    mHjacint+=(bccells[0]->dExtrapolateQuant(Varnr::mH)*=0.5*out[0]);
    mHjacint+=(bccells[1]->dExtrapolateQuant(Varnr::mH)*=-0.5*out[1]);
    mHjacint+=JacCol(bccells[0]->mHbc(),-out[0]*eye);
    jac+=mHjacint;
    eqnr+=Ns;

    // Fourth equation: continuity of total enthalpy
    JacRow conHjac(eqnr,5);
    conHjac+=JacCol(bccells[0]->Tbc(),cp*eye);
    conHjac+=JacCol(bccells[1]->Tbc(),-cp*eye);
    conHjac+=JacCol(bccells[0]->ubc(),fDFT*diagmat(u0tbc)*iDFT);    
    conHjac+=JacCol(bccells[1]->ubc(),-fDFT*diagmat(u1tbc)*iDFT);    
    jac+=conHjac;
    eqnr+=Ns;

    // Fifth equation: heat out of one segment goes into the other
    JacRow Qintjac(eqnr,5);
    Qintjac+=(bccells[0]->dExtrapolateQuant(Varnr::Q)*=out[0]);
    Qintjac+=(bccells[1]->dExtrapolateQuant(Varnr::Q)*=out[1]);
    jac+=Qintjac;
    eqnr+=Ns;

    // 6th equation: Change in exergy due to minor loss
    TRACE(15," 6th equation: Change in exergy due to minor loss");

    JacRow Exjac(eqnr,6);
    d one_eight=1.0/8.0;
    // Exjac+=JacCol(bccells[1]->pbc(),eye);
    // Exjac+=JacCol(bccells[0]->pbc(),-eye);
    Exjac+=JacCol(bccells[0]->pbc(),fDFT*diagmat(1/p0tbc)*iDFT);
    Exjac+=JacCol(bccells[1]->pbc(),-fDFT*diagmat(1/p1tbc)*iDFT);
    Exjac+=JacCol(bccells[1]->Tbc(),fDFT*diagmat(gamma/(T1tbc*(gamma-1)))*iDFT);
    Exjac+=JacCol(bccells[0]->Tbc(),-fDFT*diagmat(gamma/(T0tbc*(gamma-1)))*iDFT);
    jac+=Exjac;
  }
  void SimpleTubeConnector::setEqNrs(us firsteqnr){
  TRACE(15,"SimpleTubeConnector::setEqNrs()");
  this->firsteqnr=firsteqnr;
}
us SimpleTubeConnector::getNEqs() const {
  TRACE(15,"SimpleTubeConnector::getNEqs()");
  return Neq*Ns;
}
void SimpleTubeConnector::show(us) const{
  TRACE(15,"SimpleTubeConnector::show()");
  checkInit();
    
  cout << "SimpleTubeConnector which connects tube " << segnrs[0] <<
    " at the " << posWord(pos[0]) << " side to tube " << segnrs[1] <<
    " on the " << posWord(pos[1]) << " side.\n";
}
void SimpleTubeConnector::updateNf(){
  TRACE(15,"SimpleTubeConnector::updateNf()");
}

} // namespace tube
//////////////////////////////////////////////////////////////////////



