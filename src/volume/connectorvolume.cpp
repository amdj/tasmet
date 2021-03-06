// connectorvolume.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "connectorvolume.h"
#include "jacobian.h"
#include "staticmsg.h"
#include "tasystem.h"
#include "piston.h"
#include "duct.h"
#include "bccell.h"

#define fDFT (gc->fDFT())
#define iDFT (gc->iDFT())
#define DDTfd (gc->DDTfd())

#define Ns (gc->Ns())
#define Nf (gc->Nf())
#define eye (eye(Ns,Ns))


namespace duct {
  SPOILNAMESPACE
  common::StaticMsg<> msg;

  typedef tasystem::Jacobian jac;
  using tasystem::TaSystem;
  using tasystem::JacRow;
  using tasystem::JacCol;
  using tasystem::Jacobian;
  using tasystem::var;
  using tasystem::Globalconf;

  void DuctConnection::setPtr(const TaSystem& sys){
    TRACE(15,"DuctConnection::getPtr()");
    t=&duct::asDuct(*sys.getSeg(segid));
    c=&(t->bcCell(position));
  }

  ConnectorVolume::ConnectorVolume(d volume,bool arbitrateMass):
    volume(volume),
    arbitrateMass(arbitrateMass)
  {
    TRACE(15,"ConnectorVolume::ConnectorVolume()");
    if(volume<=0)
      throw MyError(msg("Given volume invalid, should be larger than 0. Given volume: %0.2f",volume));
  }
  ConnectorVolume::~ConnectorVolume(){
    TRACE(0,"ConnectorVolume::~ConnectorVolume()");
  }
  ConnectorVolume::ConnectorVolume(const tasystem::TaSystem& sys,\
				   const ConnectorVolume& other):
    Seg(other,sys),
    volume(other.volume),
    ductConnections(other.ductConnections),
    pistonConnections(other.pistonConnections),
    arbitrateMass(other.arbitrateMass)
  {
    TRACE(15,"ConnectorVolume::ConnectorVolume(copy)");
    // Initialize all variables
    p_=var(*gc);
    T_=var(*gc,gc->T0());
    rho_=var(*gc,gc->rho0());

    // Check if list of connections is not empty
    if(ductConnections.size()+pistonConnections.size()==0)
      throw MyError("No connections made with this volume. At least one connection is required");
    
    // Check if all pistons are really pistons
    // for(Connection& con:pistonConnections){
    //   try{
    // 	const mech::Piston&p=sys.getPiston(con.segid);
    // 	p.setConnected(con.position);
    //   }
    //   catch(...){
    // 	throw MyError(msg("Error: segment %d is not of type Piston!",con.segid));
    //   }
    // }

    // Check if all ducts are really ducts
    for(DuctConnection& con:ductConnections){
      try{
	con.setPtr(sys);
      }
      catch(...){
	throw MyError(msg("Error: segment %s is not of type Duct!",con.segid.c_str()));
      }
      
    }
    
  }
  void ConnectorVolume::addDuct(const string& segid,Pos pos){
    TRACE(15,"ConnectorVolume::addDuct()");
    ductConnections.push_back(DuctConnection(segid,pos));
  }
  void ConnectorVolume::addPiston(const string& segid,Pos pos){
    TRACE(15,"ConnectorVolume::addPiston()");
    WARN("Not yet working!");
    assert(false);
    // pistonConnections.push_back(Connection(segid,pos));    
  }
  int ConnectorVolume::arbitrateMassEq() const {
    if(arbitrateMass)
      return firsteqnr;
    else
      return -1;
  }
  void ConnectorVolume::setEqNrs(us firsteqnr){
    TRACE(50,"us ConnectorVolume::setEqNrs()");
    this->firsteqnr=firsteqnr;
  }
  void ConnectorVolume::setDofNrs(us firstdof) {
    TRACE(15,"void ConnectorVolume::setDofNrs()");
    us dof=firstdof;

    rho_.setDofNr(dof); 
    T_.setDofNr(dof+Ns); 
    p_.setDofNr(dof+2*Ns);
    
  }
  us ConnectorVolume::getNEqs() const {
    TRACE(15,"ConnectorVolume::getNEqs()");

    // conservation of mass and energy for
    // the volume. And perfectgas equation
    // of state
    us neqs=3*Ns;
    // 3 equations for each duct connected to this volume
    neqs+=3*Ns*ductConnections.size();

    return neqs;
  }
  us ConnectorVolume::getNDofs() const {
    TRACE(15,"us ConnectorVolume::getNDofs()");
    return 3*Ns;
  }
  void ConnectorVolume::setRes(const vd& res) {
    TRACE(15,"void ConnectorVolume::setRes()");
    assert(res.size()==getNDofs());
    us dofnr=0;
    rho_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    T_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    p_.setadata(res.subvec(dofnr,dofnr+Ns-1)); 
  }
  vd ConnectorVolume::getRes() const {
    TRACE(15,"vd ConnectorVolume::getRes()");
    vd res(getNDofs());
    us dofnr=0;
    res.subvec(dofnr,dofnr+Ns-1)=rho_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=T_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=p_(); 
    return res;
  }

  // Definitions
  #define Definitions				\
    const d p0=gc->p0();			\
    const d T0=gc->T0();			\
    const d cp=gc->gas().cp(T0);		\
    const d rho0=gc->gas().rho(T0,p0);		\
    const d Rs=gc->gas().Rs();			\
    const d gamma=gc->gas().gamma(T0);		\
    const d gamfac=1/(gamma-1);			\
    const d gampow=gamma/(gamma-1);		\
    const vd& rhot=rho_.tdata();		\
    const vd pt=p_.tdata()+p0;			\
    const vd& Tt=T_.tdata();

  #define DuctDefinitions				\
    const vd& utbc=con.c->ubc().tdata();		\
    const vd& Ttbc=con.c->Tbc().tdata();		\
    vd half_ubcsqt_div_cpT=pow(utbc,2)/(2*cp*Ttbc);	\
    const vd ptbc=con.c->pbc().tdata()+p0;		\
    const  vd facM0=pow(1+half_ubcsqt_div_cpT,gampow);			\
    const vd diff_facM0=gampow*pow(1+half_ubcsqt_div_cpT,gampow-1);	

  vd ConnectorVolume::error() const {
    TRACE(15,"void ConnectorVolume::error()");
    
    vd error(getNEqs());
    us eqnr=0;

    Definitions

    // ************************************************************
    // Mass conservation equation
    error.subvec(eqnr,eqnr+Ns-1)=DDTfd*(rho_()*volume);
    // Mass flow contribution from connected ducts
    for(const DuctConnection& con: ductConnections){
      assert(con.c&&con.t);
      // Flow out of the volume is defined positive
      d sign=(con.position==Pos::left?1:-1);
      error.subvec(eqnr,eqnr+Ns-1)+=sign*con.c->mbc()();
    }

    eqnr+=Ns;
    
    // ************************************************************
    // Energy conservation equation
    error.subvec(eqnr,eqnr+Ns-1)=DDTfd*(gamfac*p_()*volume);
    // Energy flow contribution from connected ducts
    for(const DuctConnection& con: ductConnections){
      assert(con.c&&con.t);
      // Flow out of the volume is defined positive, hence flow out of
      // tube should be subtracted
      d sign=(con.position==Pos::left?-1:1);
      error.subvec(eqnr,eqnr+Ns-1)+=-sign*con.c->mHbc()();
      error.subvec(eqnr,eqnr+Ns-1)+=-sign*con.c->extrapolateQuant(Varnr::Q);
    }

    eqnr+=Ns;    

    // ************************************************************
    // Equation of state
    error.subvec(eqnr,eqnr+Ns-1)=p_()-fDFT*(Tt%rhot*Rs);
    error(eqnr)+=p0;

    eqnr+=Ns;

    // ************************************************************
    // Duct equations
    for(const DuctConnection& con: ductConnections){
      assert(con.c&&con.t);
      DuctDefinitions
      // First duct equation: extrapolated mH equals mH
      error.subvec(eqnr,eqnr+Ns-1)=con.c->extrapolateQuant(Varnr::mH)
	-con.c->mHbc()();
      eqnr+=Ns;

      // Second equation: total pressure in duct equals total pressure
      // in volume
      error.subvec(eqnr,eqnr+Ns-1)=fDFT*(ptbc%facM0-pt);
      eqnr+=Ns;

      // Third equation: conduction from duct to volume
      // Compute conduction path length
      d Lcon=pow(volume,1.0/3.0);
      // Compute cross sectional area of duct exit
      d Sf=con.c->Sfbc();
      // Kappa of volume in time domain.
      vd kappat=gc->gas().kappa(Tt);

      // Watch it! This sign is DIFFERENT from the one above!!! Here
      // OUTFLOW from the DUCT is defined POSITIVE
      d sign=(con.position==Pos::left?-1:1);

      error.subvec(eqnr,eqnr+Ns-1)=\
	// Heat flow from volume to duct
	fDFT*((Sf/Lcon)*kappat%(Tt-con.c->Tbc().tdata()))+	\
	// Heat flow from the duct to the volume
	sign*con.c->extrapolateQuant(Varnr::Q);

      eqnr+=Ns;
    } // Duct equations

    return error;
  }
  void ConnectorVolume::jac(Jacobian& jac) const{
    TRACE(15,"void ConnectorVolume::jac()");

    us eqnr=firsteqnr;

    Definitions

    // ************************************************************
    // Mass conservation equation
    JacRow continuityeq(eqnr,1+ductConnections.size()+pistonConnections.size());
    continuityeq+=JacCol(rho_,volume*DDTfd);
    // Mass flow contribution from connected ducts
    for(const DuctConnection& con: ductConnections){
      assert(con.c&&con.t);
      // Flow out of the volume is defined positive
      d sign=(con.position==Pos::left?1:-1);
      continuityeq+=JacCol(con.c->mbc(),sign*eye);
    }
    // Mass flow contribution from connected pistons

    jac+=continuityeq;
    eqnr+=Ns;

    // ************************************************************
    // Energy conservation equation
    JacRow energyeq(eqnr,1+ductConnections.size()+pistonConnections.size());
    energyeq+=JacCol(p_,volume*gamfac*DDTfd);
    // Energy flow contribution from connected ducts
    for(const DuctConnection& con: ductConnections){
      assert(con.c&&con.t);
      // Flow out of the volume is defined positive
      d sign=(con.position==Pos::left?1:-1);
      energyeq+=JacCol(con.c->mHbc(),sign*eye);
      energyeq+=(con.c->dExtrapolateQuant(Varnr::Q)*=sign);
    }
    // Energy flow contribution from connected pistons

    jac+=energyeq;
    eqnr+=Ns;

    // ************************************************************
    // Equation of state 
    JacRow eos(eqnr,3);
    eos+=JacCol(p_,eye);
    eos+=JacCol(rho_,-fDFT*diagmat(Tt*Rs)*iDFT);
    eos+=JacCol(T_,-fDFT*diagmat(rhot*Rs)*iDFT);
    jac+=eos;
    eqnr+=Ns;

    // ************************************************************
    // Duct equations
    for(const DuctConnection& con: ductConnections){
      assert(con.c&&con.t);
      DuctDefinitions
      // First duct equation: extrapolated mH equals mH
      JacRow mHismH(eqnr,5);
      mHismH+=con.c->dExtrapolateQuant(Varnr::mH);
      mHismH+=JacCol(con.c->mHbc(),-eye);
      jac+=mHismH;
      eqnr+=Ns;

      // Second equation: total pressure in duct equals total pressure
      // in volume
      JacRow preseq(eqnr,2);
      preseq+=JacCol(con.c->pbc(),fDFT*diagmat(facM0)*iDFT);
      preseq+=JacCol(con.c->ubc(),fDFT*\
		     diagmat(diff_facM0%(utbc/(cp*Ttbc)))	\
		     *iDFT);      
      preseq+=JacCol(p_,-eye);
      jac+=preseq;
      eqnr+=Ns;

      // Third equation: heat flow out of Duct equals heat flow from
      // volume to duct

      // Compute conduction path length
      d Lcon=pow(volume,1.0/3.0);
      // Compute cross sectional area of duct exit
      d Sf=con.c->Sfbc();
      // Kappa of volume in time domain.
      vd kappat=gc->gas().kappa(Tt);

      // Watch it! This sign is DIFFERENT from the one above!!! Here
      // OUTFLOW from the DUCT is defined POSITIVE
      d sign=(con.position==Pos::left?-1:1);

      JacRow condeq(eqnr,4);
      condeq+=JacCol(T_,fDFT*(Sf/Lcon)*diagmat(kappat)*iDFT);
      condeq+=JacCol(con.c->Tbc(),-fDFT*(Sf/Lcon)*diagmat(kappat)*iDFT);
      condeq+=(con.c->dExtrapolateQuant(Varnr::Q)*=sign);

      jac+=condeq;
      eqnr+=Ns;

    } // Duct equations
  }
  void ConnectorVolume::show(us detailnr) const {
    TRACE(15,"void ConnectorVolume::show()");
    cout << "ConnectorVolume segment with volume: "<< volume<<"\n";
  }
  d ConnectorVolume::getMass() const {
    TRACE(15,"void ConnectorVolume::getMass()");
    return rho_(0)*volume;
  }
  void ConnectorVolume::dmtotdx(vd& dmtotdx) const {
    TRACE(15,"void ConnectorVolume::dmtotdx()");
    us rhodof=rho_.getDofNr();
    dmtotdx(rhodof)=volume;
  }
  void ConnectorVolume::domg(vd& domg) const{
    TRACE(15,"void ConnectorVolume::domg()");
    WARN("Code has not been checked");
    WARN("Needs implementation");    
    assert(false);
  }
  void ConnectorVolume::updateNf(){
    TRACE(15,"ConnectorVolume::updateNf()");
    p_.updateNf();
    T_.updateNf();
    rho_.updateNf();
  }
  void ConnectorVolume::resetHarmonics(){
    TRACE(15,"Piston::updateNf()");
    p_.resetHarmonics();
    T_.resetHarmonics();
    rho_.resetHarmonics();
  }

} // namespace duct
//////////////////////////////////////////////////////////////////////
