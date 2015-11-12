// piston.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "piston.h"
#include "jacobian.h"
#include "staticmsg.h"

#define fDFT (gc->fDFT())
#define iDFT (gc->iDFT())
#define DDTfd (gc->DDTfd())

#define Ns (gc->Ns())
#define Nf (gc->Nf())
#define eye (eye(Ns,Ns))

namespace{
    common::StaticMsg<> msg;
}
namespace mech {

  SPOILNAMESPACE
  typedef tasystem::Jacobian jac;
  using tasystem::JacRow;
  using tasystem::JacCol;
  using tasystem::Jacobian;
  using tasystem::var;
  using tasystem::Globalconf;

  PistonConfiguration::PistonConfiguration(d Sl,d Sr,d V0l,d V0r,d M,d Km,d Cm,d Stl,d Str):
      M(M),Sr(Sr),Sl(Sl),Km(Km),Cm(Cm),V0l(V0l),V0r(V0r),
      Stl(Stl),Str(Str)
    {
      if(V0l<=0)
	throw MyError("Illegal value for V0l given");
      if(V0r<=0)
	throw MyError("Illegal value for V0r given");
      // If Stl is undefined, we make it a cylinder
      if(Stl<0){
        d L=V0l/Sl;
        d circumference=2*number_pi*sqrt(Sl/number_pi);
        Stl=Sl+L*circumference;
      }
      // Same for Str
      if(Str<0){
        d L=V0r/Sr;
        d circumference=2*number_pi*sqrt(Sr/number_pi);
        Str=Sr+L*circumference;
      }
    }
  Piston::Piston(const PistonConfiguration& pc,bool arbitrateMass):
      Seg(),
      pc(pc),
      arbitrateMass(arbitrateMass)
    {
      TRACE(15,"Piston()");
    }

  Piston::Piston(const tasystem::TaSystem& sys,const Piston& other):
    Seg(other,sys),
    pc(other.pc),
    T0(other.T0),
    arbitrateMass(other.arbitrateMass),
    massL(other.massL),
    massR(other.massR)
  {
    TRACE(15,"Piston::Piston()");

    // If temperature not initialized
    if(T0<=0)
      T0=gc->T0();

    // Initialize all variables
    xp_=var(*gc);
    Fp_=var(*gc);

    pl_=var(*gc);
    pr_=var(*gc);
    ml_=var(*gc);
    mr_=var(*gc);
    Tl_=var(*gc,T0);
    Tr_=var(*gc,T0);
    rhol_=var(*gc,gc->rho0());
    rhor_=var(*gc,gc->rho0());
    mHl_=var(*gc);
    mHr_=var(*gc);
    if(massL<0)
      massL=gc->rho0()*pc.V0l;
    if(massR<0)
      massR=gc->rho0()*pc.V0r;
  }
  Piston::~Piston(){}
  int Piston::arbitrateMassEq() const {
    TRACE(15,"Piston::arbitrateMassEq()");
    VARTRACE(80,firsteqnr);
    VARTRACE(80,Ns);    
    if(arbitrateMass){
      if(leftConnected)
	return firsteqnr+Ns;
      else if(rightConnected)
	return firsteqnr+4*Ns;
      else
	return -1;
    }
    else
      return -1;
  }
  void Piston::setEqNrs(us firsteqnr){
    TRACE(15,"us Piston::setEqNrs()");
    this->firsteqnr=firsteqnr;
  }
  void Piston::setDofNrs(us firstdof) {
    TRACE(15,"void Piston::setDofNrs()");
    us dof=firstdof;
    xp_.setDofNr(dof); dof+=Ns;
    Fp_.setDofNr(dof); dof+=Ns;

    rhol_.setDofNr(dof); dof+=Ns;
    Tl_.setDofNr(dof); dof+=Ns;
    pl_.setDofNr(dof); dof+=Ns;
    ml_.setDofNr(dof); dof+=Ns;
    mHl_.setDofNr(dof); dof+=Ns;
    
    rhor_.setDofNr(dof); dof+=Ns;
    Tr_.setDofNr(dof); dof+=Ns;
    pr_.setDofNr(dof); dof+=Ns;
    mr_.setDofNr(dof); dof+=Ns;
    mHr_.setDofNr(dof); dof+=Ns;
    
  }
  us Piston::getNEqs() const {
    TRACE(15,"void Piston::getNEqs()");
    // 1 for equation of motion of piston
    // 2 for mass conservation left and right volumes
    // 2 for energy conservation left and right volumes
    // 2 for equations of state
    // ------ +
    // 7 equations

    us neqs=7*Ns;
    // If not connected, one to set mass flow left and right to zero
    // and one to set enthalpy flow left and right to zero
    if(!leftConnected){
      TRACE(15,"Left side not connected");
      neqs+=2*Ns;                 // For ml=0
    }
    if(!rightConnected){
      neqs+=2*Ns;                 // For mr=0
      TRACE(15,"Right side not connected");
    }
    return neqs;
  }
  us Piston::getNDofs() const {
    TRACE(15,"us Piston::getNDofs()");
    return 12*Ns;
  }
  void Piston::setRes(const vd& res) {
    TRACE(15,"void Piston::setRes()");
    us dofnr=0;
    // order: xp,Fp,rhol,Tl,pl,ml,rhor,Tr,pr,mr
    xp_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    Fp_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;

    rhol_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    Tl_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    pl_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    ml_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    mHl_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;

    rhor_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    Tr_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    pr_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    mr_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    mHr_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
  }
  vd Piston::getRes() const {
    TRACE(15,"vd Piston::getRes()");
    vd res(getNDofs());
    us dofnr=0;
    // order: xp,Fp,rhol,Tl,pl,ml,rhor,Tr,pr,mr
    res.subvec(dofnr,dofnr+Ns-1)=xp_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=Fp_(); dofnr+=Ns;

    res.subvec(dofnr,dofnr+Ns-1)=rhol_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=Tl_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=pl_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=ml_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=mHl_(); dofnr+=Ns;

    res.subvec(dofnr,dofnr+Ns-1)=rhor_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=Tr_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=pr_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=mr_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=mHr_(); dofnr+=Ns;
    return res;
  }
  vd Piston::error() const {
    TRACE(15,"void Piston::error()");
    
    vd error(getNEqs());
    us eqnr=0;

    const d p0=gc->p0();
    const d rho0=gc->gas().rho(T0,p0);
    const d gamma=gc->gas().gamma(T0);
    const d gamfac=1/(gamma-1);


    const vd& rholt=rhol_.tdata();
    const vd& rhort=rhor_.tdata();
    const vd plt=pl_.tdata()+p0;
    const vd prt=pr_.tdata()+p0;

    const vd& Tlt=Tl_.tdata();
    const vd& Trt=Tr_.tdata();
    const vd& mlt=ml_.tdata();
    const vd& mrt=mr_.tdata();
    const vd& xpt=xp_.tdata();
    // Volume as a function of time
    const vd Vlt=pc.V0l+xpt*pc.Sl;
    const vd Vrt=pc.V0r-xpt*pc.Sr;

    // dVdt in time domain of left volume
    const vd dVldtt=iDFT*DDTfd*(xp_()*pc.Sl);
    // dVdt in time domain of right volume
    const vd dVrdtt=iDFT*DDTfd*(-xp_()*pc.Sr);

    // ************************************************************
    // Equation of motion of the piston
    error.subvec(eqnr,eqnr+Ns-1)=\
      pc.M*DDTfd*DDTfd*xp_()      // "Inertia" force
      +pc.Cm*DDTfd*xp_()          // Damping force
      +pc.Km*xp_()                // Spring force
      -pl_()*pc.Sl                // "Pushes"
      +pr_()*pc.Sr                // "Pushes back"
      -Fp_();                   // "External force on piston

    eqnr+=Ns;

    // ************************************************************
    // Mass conservation left side
    if(!leftConnected){
      assert(massL>0);
      error.subvec(eqnr,eqnr+Ns-1)=fDFT*(rholt%Vlt);
      error(eqnr)-=massL;

      // Not connected left side, so set ml to zero
      eqnr+=Ns;
      error.subvec(eqnr,eqnr+Ns-1)=ml_();
      // Not connected left side, so set mHl to zero
      eqnr+=Ns;
      error.subvec(eqnr,eqnr+Ns-1)=mHl_();
    }
    else{
      error.subvec(eqnr,eqnr+Ns-1)=\
        DDTfd*fDFT*(rholt%Vlt)+ml_();
    }
    eqnr+=Ns;

    // ************************************************************
    // Mass conservation right side
    if(!rightConnected){
      assert(massR>0);
      error.subvec(eqnr,eqnr+Ns-1)=fDFT*(rhort%Vrt);
      error(eqnr)-=massR;

      // Not connected right side, so set mr to zero
      eqnr+=Ns;
      error.subvec(eqnr,eqnr+Ns-1)=mr_();

      // Not connected right side, so set mHr to zero
      eqnr+=Ns;
      error.subvec(eqnr,eqnr+Ns-1)=mHr_();

    }
    else{
      error.subvec(eqnr,eqnr+Ns-1)=\
        DDTfd*fDFT*(rhort%Vrt)+mr_();
    }
    eqnr+=Ns;

    // ************************************************************
    // Energy consrvation left side
    error.subvec(eqnr,eqnr+Ns-1)=gamfac*DDTfd*fDFT*(plt%Vlt)
      +fDFT*(plt%dVldtt) // Work contribution to change in
                                    // energy
      +mHl_();            // Enthalpy flow out of segment

    // if(!leftConnected){
      // Overwrite time-average part with constraint on time-average
      // internal energy, by fixing the time-averaged temperature
      error(eqnr)=Tl_()(0)-T0;

      // Isentropic model
      // error.subvec(eqnr,eqnr+Ns-1)=pl_()/p0;
      // error(eqnr)+=1;
      // error.subvec(eqnr,eqnr+Ns-1)+=-fDFT*pow(rhol_.tdata()/rho0,gamma);
    // }
    eqnr+=Ns;

    // ************************************************************
    // Energy consrvation right side
    error.subvec(eqnr,eqnr+Ns-1)=gamfac*DDTfd*fDFT*(prt%Vrt)
      +fDFT*(prt%dVrdtt) // Work contribution to change in
      // energy
      +mHr_();            // Enthalpy flow out of segment

    // if(!rightConnected){
      // Overwrite time-average part with constraint on time-average
      // internal energy, by fixing the time-averaged temperature
      error(eqnr)=Tr_()(0)-T0;

      // Isentropic model
      // error.subvec(eqnr,eqnr+Ns-1)=pr_()/p0;
      // error(eqnr)+=1;
      // error.subvec(eqnr,eqnr+Ns-1)+=-fDFT*pow(rhor_.tdata()/rho0,gamma);
    // }
    eqnr+=Ns;

    // Specific gas constant
    d Rs=gc->gas().Rs();

    // ************************************************************
    // Equation of state left side
    error.subvec(eqnr,eqnr+Ns-1)=pl_()-fDFT*(Tl_.tdata()%rhol_.tdata()*Rs);
    error(eqnr)+=p0;
    eqnr+=Ns;

    // ************************************************************
    // Equation of state right side
    error.subvec(eqnr,eqnr+Ns-1)=pr_()-fDFT*(Tr_.tdata()%rhor_.tdata()*Rs);
    error(eqnr)+=p0;

    TRACE(15,"Returning from Piston::Error()");
    return error;
  }
  void Piston::jac(Jacobian& jac) const{
    TRACE(15,"void Piston::jac()");
    TRACE(15,"forsteqnr: "<< firsteqnr);

    us eqnr=firsteqnr;

    const d p0=gc->p0();
    const d rho0=gc->gas().rho(T0,p0);
    const d gamma=gc->gas().gamma(T0);
    const d gamfac=1/(gamma-1);
    const d cp=gc->gas().cp(T0);

    const vd& rholt=rhol_.tdata();
    const vd& rhort=rhor_.tdata();
    const vd plt=pl_.tdata()+p0;
    const vd prt=pr_.tdata()+p0;

    const vd& Tlt=Tl_.tdata();
    const vd& Trt=Tr_.tdata();
    const vd& mlt=ml_.tdata();
    const vd& mrt=mr_.tdata();
    const vd& xpt=xp_.tdata();
    // Volume as a function of time
    const vd Vlt=pc.V0l+xpt*pc.Sl;
    const vd Vrt=pc.V0r-xpt*pc.Sr;

    // dVdt in time domain of left volume
    const vd dVldtt=iDFT*DDTfd*(xp_()*pc.Sl);
    // dVdt in time domain of right volume
    const vd dVrdtt=iDFT*DDTfd*(-xp_()*pc.Sr);

    // ************************************************************
    // 1: Equation of motion of the piston
    JacRow eomjac(eqnr,4);
    eomjac+=JacCol(xp_,pc.M*DDTfd*DDTfd+pc.Cm*DDTfd+pc.Km*eye);
    eomjac+=JacCol(Fp_,-eye);
    eomjac+=JacCol(pl_,-pc.Sl*eye);
    eomjac+=JacCol(pr_,pc.Sr*eye);
    jac+=eomjac;

    eqnr+=Ns;

    // ************************************************************
    // 2: Mass conservation left side
    if(!leftConnected){
      JacRow mcljac(eqnr,2);
      mcljac+=JacCol(rhol_,fDFT*diagmat(Vlt)*iDFT);
      mcljac+=JacCol(xp_,fDFT*diagmat(rholt*pc.Sl)*iDFT);
      jac+=mcljac;

      // Not connected to leftside so ml to zero
      eqnr+=Ns;
      jac+=JacRow(eqnr,JacCol(ml_,eye));
      // Not connected to leftside so mHl to zero
      eqnr+=Ns;
      jac+=JacRow(eqnr,JacCol(mHl_,eye));
    }
    else{
      JacRow mcljac(eqnr,3);
      mcljac+=JacCol(rhol_,DDTfd*fDFT*diagmat(Vlt)*iDFT);
      mcljac+=JacCol(xp_,DDTfd*fDFT*diagmat(rhol_.tdata()*pc.Sl)*iDFT);
      mcljac+=JacCol(ml_,eye);
      jac+=mcljac;
    }
    eqnr+=Ns;

    // ************************************************************    
    // 3: Mass conservation right side
    if(!rightConnected){
      JacRow mcrjac(eqnr,2);
      mcrjac+=JacCol(rhor_,fDFT*diagmat(Vrt)*iDFT);
      mcrjac+=JacCol(xp_,fDFT*diagmat(-rhor_.tdata()*pc.Sr)*iDFT);
      jac+=mcrjac;

      // Not connected to right side so mr is set to zero
      eqnr+=Ns;
      jac+=JacRow(eqnr,JacCol(mr_,eye));
      // Not connected to right side so mHr is set to zero
      eqnr+=Ns;
      jac+=JacRow(eqnr,JacCol(mHr_,eye));

    }
    else{
      JacRow mcrjac(eqnr,3);
      mcrjac+=JacCol(rhor_,DDTfd*fDFT*diagmat(Vrt)*iDFT);
      mcrjac+=JacCol(xp_,DDTfd*fDFT*diagmat(-rhor_.tdata()*pc.Sr)*iDFT);
      mcrjac+=JacCol(mr_,eye);
      jac+=mcrjac;
    }
    eqnr+=Ns;

    // ************************************************************
    {
      // Energy equation left side:
      // 1/(gamma-1)*d/dt(p*Vl)+p*dVl/dt+ml*hl=0

      // 4: energy conservation left side
      // 4 columns in this row if not connected, otherwise 3
      JacRow enl(eqnr,2+(!leftConnected?1:0));
    
      // Terms due to change in internal energy
      dmat enlmat_pl=gamfac*DDTfd*fDFT*diagmat(Vlt)*iDFT;
      dmat enlmat_x=gamfac*DDTfd*fDFT*diagmat(plt*pc.Sl)*iDFT;

      // Terms due to work done on fluid
      enlmat_pl+=fDFT*diagmat(dVldtt)*iDFT;
      enlmat_x+=fDFT*diagmat(pc.Sl*plt)*iDFT*DDTfd;
    
      // if(!leftConnected) {
        // Overwrite time-average part with constraint on time-average
        // internal energy, by fixing the time-averaged temperature
        dmat enlmat_T(Ns,Ns,fillwith::zeros);
        enlmat_T(0,0)=1;

        // Zero out first row
        enlmat_pl.row(0).zeros();
        enlmat_x.row(0).zeros();

        // Add Jacobian terms corresponding to left temperature
        enl+=JacCol(Tl_,enlmat_T);          

      // }

      enl+=JacCol(mHl_,eye);
      enl+=JacCol(pl_,enlmat_pl);
      enl+=JacCol(xp_,enlmat_x);    
      jac+=enl;
      eqnr+=Ns;
    }
    
    // ************************************************************
    {
      // Energy equation right side
      // 1/(gamma-1)*d/dt(p*Vr)+p*dVr/dt+mr*hr=0
      // 3 columns in this row if not connected, otherwise 2
      JacRow enr(eqnr,2+(!leftConnected?1:0));
    
      // Terms due to change in internal energy
      dmat enrmat_pr=gamfac*DDTfd*fDFT*diagmat(Vrt)*iDFT;
      dmat enrmat_x=gamfac*DDTfd*fDFT*diagmat(-prt*pc.Sr)*iDFT;

      // Terms due to work done on fluid
      enrmat_pr+=fDFT*diagmat(dVrdtt)*iDFT;
      enrmat_x+=-fDFT*diagmat(pc.Sr*prt)*iDFT*DDTfd;
    
      // if(!rightConnected) {
        // Overwrite time-average part with constraint on time-average
        // internal energy, by fixing the time-averaged temperature
        dmat enrmat_T(Ns,Ns,fillwith::zeros);
        enrmat_T(0,0)=1;

        // Zero out first row
        enrmat_pr.row(0).zeros();
        enrmat_x.row(0).zeros();


        // Add Jacobian terms corresponding to left temperature
        enr+=JacCol(Tr_,enrmat_T);          

      // }

      enr+=JacCol(mHr_,eye);
      enr+=JacCol(pr_,enrmat_pr);
      enr+=JacCol(xp_,enrmat_x);    
      jac+=enr;
      eqnr+=Ns;
    }
    
    // ************************************************************
    const d Rs=gc->gas().Rs();
    // Equation of state left side
    {
      JacRow eosl(eqnr,3);
      eosl+=JacCol(pl(),eye);
      eosl+=JacCol(rhol_,-fDFT*diagmat(Tl_.tdata()*Rs)*iDFT);
      eosl+=JacCol(Tl_,-fDFT*diagmat(rhol_.tdata()*Rs)*iDFT);
      jac+=eosl;
      eqnr+=Ns;
    }
    // ************************************************************
    // Equation of state right side
    {
      JacRow eosr(eqnr,3);
      eosr+=JacCol(pr(),eye);
      eosr+=JacCol(rhor_,-fDFT*diagmat(Tr_.tdata()*Rs)*iDFT);
      eosr+=JacCol(Tr_,-fDFT*diagmat(rhor_.tdata()*Rs)*iDFT);
      jac+=eosr;
      eqnr+=Ns;
    }
  }
  var Piston::upiston() const{
    return var(*gc,DDTfd*xpiston()());
  }
  void Piston::show(us detailnr) const {
    TRACE(15,"void Piston::show()");
    cout << "Piston segment:\n"
         << "Piston mass: " << pc.M
         << "\nPiston mechanical stiffness: "<< pc.Km
         << "\nPiston mechanical damping: "<< pc.Cm
         << "\nLeft volume: "<< pc.V0l    
         << "\nRight volume: "<< pc.V0r << endl;    
    if(!leftConnected)
      cout << "Left initial mass: "<< massL << endl;    
    if(!rightConnected)
      cout << "Right initial mass: "<< massR << endl;    
  }
  d Piston::getMass() const {
    TRACE(15,"void Piston::getMass()");
    // Only the mass counts which is connected to other segments
    d mass=0;
    if(leftConnected){
      mass+=arma::dot( fDFT.row(0) , rhol_.tdata()%(pc.V0l+pc.Sl*xp_.tdata()) );
    }
    if(rightConnected){
      mass+=arma::dot( fDFT.row(0) , rhor_.tdata()%(pc.V0r-pc.Sr*xp_.tdata()) );
    }
    return mass;
  }
  void Piston::dmtotdx(vd& dmtotdx) const {
    TRACE(15,"void Piston::dmtotdx()");
    if(leftConnected){
      us rhodof=rhol_.getDofNr();
      us xdof=xp_.getDofNr();
      dmtotdx.subvec(rhodof,rhodof+Ns-1)=(fDFT.row(0)*			\
					  diagmat(pc.V0l+pc.Sl*xp_.tdata()*iDFT)).t();
      dmtotdx.subvec(xdof,xdof+Ns-1)=(fDFT.row(0)*			\
					  diagmat(pc.Sl*rhol_.tdata()*iDFT)).t();
    }
    if(rightConnected){
      us rhodof=rhor_.getDofNr();
      us xdof=xp_.getDofNr();
      dmtotdx.subvec(rhodof,rhodof+Ns-1)=(fDFT.row(0)*			\
					  diagmat(pc.V0r-pc.Sr*xp_.tdata())*iDFT).t();
      dmtotdx.subvec(xdof,xdof+Ns-1)=(fDFT.row(0)*			\
				      diagmat(-pc.Sr*rhor_.tdata())*iDFT).t();
    }

  }
  void Piston::domg(vd& domg) const{
    TRACE(15,"void Piston::domg()");
    WARN("Code has not been checked");
    WARN("Needs implementation");    
    assert(false);
  }
  void Piston::updateNf(){
    TRACE(15,"Piston::updateNf()");
    xp_.updateNf();
    Fp_.updateNf();    
    pl_.updateNf();
    pr_.updateNf();
    ml_.updateNf();
    mr_.updateNf();
    Tl_.updateNf();
    Tr_.updateNf();
    rhol_.updateNf();
    rhor_.updateNf();
  }
  void Piston::resetHarmonics(){
    TRACE(15,"Piston::updateNf()");
    xp_.resetHarmonics();
    Fp_.resetHarmonics();
    pl_.resetHarmonics();
    pr_.resetHarmonics();
    ml_.resetHarmonics();
    mr_.resetHarmonics();
    Tl_.resetHarmonics();
    Tr_.resetHarmonics();
    rhol_.resetHarmonics();
    rhor_.resetHarmonics();
  }
  var Piston::Vl() const {
    TRACE(15,"var Piston::Vl()");
    return var(*gc,pc.V0l+xp_.tdata()*pc.Sl,false);
  }
  var Piston::Vr() const {
    TRACE(15,"var Piston::Vr()");
    return var(*gc,pc.V0r-xp_.tdata()*pc.Sr,false);
  }
  bool Piston::isConnected(Pos side) const {
    TRACE(15,"bool Piston::isConnected()");
    return side==Pos::left?leftConnected:rightConnected;
  }
  void Piston::setConnected(Pos side) const {
    TRACE(15,"void Piston::setConnected()");
    bool& connected=(side==Pos::left?leftConnected:rightConnected);
    if(connected)
      throw MyError(msg("Error: piston is already connected on side %s.",duct::posWord(side)));
    connected=true;
  }
} // namespace mech
//////////////////////////////////////////////////////////////////////
