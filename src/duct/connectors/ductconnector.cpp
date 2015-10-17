// ductconnector.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
// Connect two duct-type segments by using conservation of mass,
// momentum and energy. Number of equations 
//////////////////////////////////////////////////////////////////////

#include "ductconnector.h"
#include "tasystem.h"
#include "duct.h"
#include "bccell.h"
#include "constants.h"
#include "jacobian.h"
#include "energy.h"

#define fDFT (gc->fDFT())
#define iDFT (gc->iDFT())
#define DDTfd (gc->DDTfd())

#define Ns (gc->Ns())
#define eye (arma::eye(Ns,Ns))

namespace duct {
  using tasystem::TaSystem;
  using tasystem::JacRow;
  using tasystem::JacCol;

  DuctConnector::DuctConnector(const string& seg1,Pos pos1,\
			       const string& seg2,Pos pos2,d K1to2,d K2to1):
    K1to2(K1to2),
    K2to1(K2to1)
  {
    TRACE(15,"DuctConnector::DuctConnector()");
  
    segids[0]=seg1;
    segids[1]=seg2;
    pos[0]=pos1;
    pos[1]=pos2;
      
  }
  DuctConnector::DuctConnector(const DuctConnector& o,
			       const TaSystem& sys):
    Connector(o,sys),
    segids(o.segids),
    pos(o.pos),
    K1to2(o.K1to2),
    K2to1(o.K2to1)
  {
    
    bccells[0]=&sys.getDuct(segids[0]).bcCell(pos[0]);
    bccells[1]=&sys.getDuct(segids[1]).bcCell(pos[1]);
    assert(bccells[0]&&bccells[1]);

    if(pos[0]==Pos::left) 
      out[0]=-1;
    if(pos[1]==Pos::left)
      out[1]=-1;
    
  }

  #define Defs							\
    d T0=gc->T0();						\
    d gamma=gc->gas().gamma(T0);				\
    d p0=gc->p0();						\
    d cp=gc->gas().cp(T0);					\
    d Rs=gc->gas().Rs();					\
    d Sf0=bccells[0]->Sfbc();					\
    d Sf1=bccells[1]->Sfbc();					\
    const vd p0tbc=bccells[0]->pbc().tdata()+p0;		\
    const vd p1tbc=bccells[1]->pbc().tdata()+p0;		\
    const vd& T0tbc=bccells[0]->Tbc().tdata();			\
    const vd& T1tbc=bccells[1]->Tbc().tdata();			\
    const vd& m0tbc=bccells[0]->mbc().tdata();			\
    const vd& m1tbc=bccells[1]->mbc().tdata();			\
    const vd& T0t=bccells[0]->T().tdata();			\
    const vd& T1t=bccells[1]->T().tdata();			\
    const vd& rho0tbc=bccells[0]->rhobc().tdata();		\
    const vd& rho1tbc=bccells[1]->rhobc().tdata();		\
    const vd& u0tbc=bccells[0]->ubc().tdata();			\
    const vd& u1tbc=bccells[1]->ubc().tdata();			\
    vd half_u0sqt_div_cpT=pow(u0tbc,2)/(2*cp*T0tbc);		\
    vd half_u1sqt_div_cpT=pow(u1tbc,2)/(2*cp*T1tbc);		\
    d gampow=gamma/(gamma-1);					\
    vd facM0=pow(1+half_u0sqt_div_cpT,gampow);			\
    vd facM1=pow(1+half_u1sqt_div_cpT,gampow);			\
    vd diff_facM0=gampow*pow(1+half_u0sqt_div_cpT,gampow-1);	\
    vd diff_facM1=gampow*pow(1+half_u1sqt_div_cpT,gampow-1);			

  vd DuctConnector::error() const{
    TRACE(10,"DuctConnector::error()");

    us nr=0;
    vd error(getNEqs());
    // VARTRACE(60,getNEqs());
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
    // enthalpy flows of both ducts
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


    // 5th equation: Change in total pressure over density due to minor loss
    d one_eight=1.0/8.0;

    vd ptot1=fDFT*(p1tbc%facM1);
    vd ptot0=fDFT*(p0tbc%facM0);
    error.subvec(Ns*nr,Ns*(nr+1)-1)=\
      (ptot1-ptot0)/p0;
    // vd minus_minorLoss=K1to2*(T0/T1t)*one_eight	\
    //   *pow(out[0]*abs(u1t)+u1t,2)-\
    //   K2to1*(T0/T2t)*one_eight\
    //   *pow(out[1]*abs(u2t)+u2t,2);
    // vd minus_minorLoss=zeros(Ns);
    nr++;

    // 6th eq: continuity of total heat flow
    error.subvec(Ns*nr,Ns*(nr+1)-1)=		\
      out[0]*bccells[0]->extrapolateQuant(Varnr::Qs)+	\
      out[0]*bccells[0]->extrapolateQuant(Varnr::Q)+	\
      out[1]*bccells[1]->extrapolateQuant(Varnr::Q)+	\
      out[1]*bccells[1]->extrapolateQuant(Varnr::Qs);
    nr++;

    // Solid equations ***********************************************
    bool sol0=bccells[0]->getDuct().hasSolid();
    bool sol1=bccells[1]->getDuct().hasSolid();
    if(sol0&&sol1) {
      // Simple case: continuity of temperature and solid heat flow
      error.subvec(Ns*nr,Ns*(nr+1)-1)=			\
	bccells[0]->Tsbc()()-bccells[1]->Tsbc()();
      nr++;
      error.subvec(Ns*nr,Ns*(nr+1)-1)=			\
	out[0]*bccells[0]->extrapolateQuant(Varnr::Qs)+\
	out[1]*bccells[1]->extrapolateQuant(Varnr::Qs);
      return error;
    }
    else if(sol0) { // Only Duct 0 has solid

      d Ss=bccells[0]->vSs;
      d rhs=bccells[0]->vrh*Ss/bccells[0]->vSf;
      d kappaf=gc->gas().kappa(bccells[1]->Tbc()(0));
      d h=Nu2*kappaf/rhs;
      error.subvec(Ns*nr,Ns*(nr+1)-1)=		\
	h*Ss*(bccells[0]->Tsbc()()-bccells[1]->Tbc()())		\
	-out[0]*bccells[0]->extrapolateQuant(Varnr::Qs);

      return error;
    } // Only Duct 0 has solid
    else if(sol1) { // Only Duct 1 has solid

      d Ss=bccells[1]->vSs;
      d rhs=bccells[1]->vrh*Ss/bccells[1]->vSf;
      d kappaf=gc->gas().kappa(bccells[0]->Tbc()(0));
      d h=Nu1*kappaf/rhs;
      error.subvec(Ns*nr,Ns*(nr+1)-1)=		\
	h*Ss*(bccells[1]->Tsbc()()-bccells[1]->Tbc()())		\
	-out[1]*bccells[1]->extrapolateQuant(Varnr::Qs);

      return error;
    }      // Only Duct 1 has a solid
    return error;
  }
  void DuctConnector::jac(tasystem::Jacobian& jac) const {
    TRACE(15,"Ductconnector::jac()");
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
    // // enthalpy flows of both ducts
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


    // 5th equation: Change in total pressure over density due to
    // minor loss

    JacRow Minor(eqnr,8);
    d one_eight=1.0/8.0;

    Minor+=JacCol(bccells[0]->pbc(),fDFT*diagmat(-(1/p0)*facM0)*iDFT);
    Minor+=JacCol(bccells[0]->ubc(),fDFT*diagmat(-(1/p0)*p0tbc%		\
						 diff_facM0%(u0tbc/(cp*T0tbc)))*iDFT);
    Minor+=JacCol(bccells[0]->Tbc(),fDFT*diagmat((1/p0)*p0tbc%		\
						 diff_facM0%(0.5*pow(u0tbc,2)/(cp*pow(T0tbc,2))))*iDFT);

    Minor+=JacCol(bccells[1]->pbc(),fDFT*diagmat((1/p0)*facM1)*iDFT);
    Minor+=JacCol(bccells[1]->ubc(),fDFT*diagmat((1/p0)*p1tbc%		\
						 diff_facM1%(u1tbc/(cp*T1tbc)))*iDFT);
    Minor+=JacCol(bccells[1]->Tbc(),fDFT*diagmat(-(1/p0)*p1tbc%		\
						 diff_facM1%(0.5*pow(u1tbc,2)/(cp*pow(T1tbc,2))))*iDFT);

    jac+=Minor;
    eqnr+=Ns;

    // 6th eq: continuity of total heat flow
    JacRow Qintjac(eqnr,5);
    Qintjac+=(bccells[0]->dExtrapolateQuant(Varnr::Q)*=out[0]);
    Qintjac+=(bccells[0]->dExtrapolateQuant(Varnr::Qs)*=out[0]);
    Qintjac+=(bccells[1]->dExtrapolateQuant(Varnr::Q)*=out[1]);
    Qintjac+=(bccells[1]->dExtrapolateQuant(Varnr::Qs)*=out[1]);
    jac+=Qintjac;
    eqnr+=Ns;

    // Solid equations ***********************************************
    bool sol0=bccells[0]->getDuct().hasSolid();
    bool sol1=bccells[1]->getDuct().hasSolid();
    if(sol0&&sol1) {
      // Simple case: continuity of temperature and solid heat flow
      JacRow Ts_is_Ts(eqnr,2);
      Ts_is_Ts+=JacCol(bccells[0]->Tsbc(),eye);
      Ts_is_Ts+=JacCol(bccells[1]->Tsbc(),-eye);	
      jac+=Ts_is_Ts;
      eqnr+=Ns;

      JacRow Qs_is_Qs(eqnr,4);
      Qs_is_Qs+=(bccells[0]->dExtrapolateQuant(Varnr::Qs)*=out[0]);
      Qs_is_Qs+=(bccells[1]->dExtrapolateQuant(Varnr::Qs)*=out[1]);
      jac+=Qs_is_Qs;
      return;
    }
    else if(sol0) { // Only Duct 0 has solid
      d Ss=bccells[0]->vSs;
      d rhs=bccells[0]->vrh*Ss/bccells[0]->vSf;
      d kappaf=gc->gas().kappa(bccells[1]->Tbc()(0));
      d h=Nu2*kappaf/rhs;
      JacRow Qseq(eqnr,4);
      Qseq+=JacCol(bccells[0]->Tsbc(),h*Ss*eye);
      Qseq+=JacCol(bccells[1]->Tbc(),-h*Ss*eye);
      Qseq+=(bccells[0]->dExtrapolateQuant(Varnr::Qs)*=-out[0]);
      jac+=Qseq;
      // eqnr+=Ns;			// Not required anymore
      return;
    } // Only Duct 0 has solid
    else if(sol1) { // Only Duct 1 has solid
      d Ss=bccells[1]->vSs;
      d rhs=bccells[1]->vrh*Ss/bccells[1]->vSf;
      d kappaf=gc->gas().kappa(bccells[0]->Tbc()(0));
      d h=Nu1*kappaf/rhs;
      JacRow Qseq(eqnr,4);
      Qseq+=JacCol(bccells[1]->Tsbc(),h*Ss*eye);
      Qseq+=JacCol(bccells[0]->Tbc(),-h*Ss*eye);
      Qseq+=(bccells[1]->dExtrapolateQuant(Varnr::Qs)*=-out[1]);
      jac+=Qseq;
      // eqnr+=Ns;			// Not required anymore
    } // Only Duct 1 has solid
  }
  void DuctConnector::setEqNrs(us firsteqnr){
    TRACE(15,"DuctConnector::setEqNrs()");
    this->firsteqnr=firsteqnr;
  }
  us DuctConnector::getNEqs() const {
    TRACE(15,"DuctConnector::getNEqs()");
    us neq=6;
    for(us i=0;i<2;i++)
      if(bccells[i]->getDuct().hasSolid())
	neq++;
    return neq*Ns;
  }
  void DuctConnector::show(us) const{
    TRACE(15,"DuctConnector::show()");
    
    cout << "DuctConnector which connects duct " << segids[0] <<
      " at the " << posWord(pos[0]) << " side to duct " << segids[1] <<
      " on the " << posWord(pos[1]) << " side.\n";
  }
  void DuctConnector::updateNf(){
    TRACE(15,"DuctConnector::updateNf()");
  }

} // namespace duct
//////////////////////////////////////////////////////////////////////



