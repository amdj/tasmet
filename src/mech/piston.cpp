// piston.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "piston.h"
#include "jacobian.h"

#define fDFT (gc->fDFT)
#define iDFT (gc->iDFT)
#define DDTfd (gc->DDTfd)

#define Ns (gc->Ns())
#define Nf (gc->Nf())
#define eye (eye(Ns,Ns))


namespace mech {
  SPOILNAMESPACE
  typedef tasystem::Jacobian jac;
  using tasystem::JacRow;
  using tasystem::JacCol;
  using tasystem::Jacobian;
  using tasystem::var;
  using tasystem::Globalconf;

  Piston::Piston(const tasystem::TaSystem& sys,const Piston& other):
    Seg(other,sys),
    M(other.M),
    Sr(other.Sr),
    Sl(other.Sl),
    Km(other.Km),
    Cm(other.Cm),
    V0l(other.V0l),
    V0r(other.V0r),
    Stl(other.Stl),
    Str(other.Str),
    T0(other.T0)
  {
    TRACE(15,"Piston::Piston()");

    // If temperature not initialized
    if(T0<=0)
      T0=gc->T0();

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

    
    if(!leftConnected && massL<0)
      massL=rhol_(0)*V0l;
    if(!rightConnected && massR<0)
      massR=rhor_(0)*V0r;
  }
  Piston::~Piston(){}
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

    rhor_.setDofNr(dof); dof+=Ns;
    Tr_.setDofNr(dof); dof+=Ns;
    pr_.setDofNr(dof); dof+=Ns;
    mr_.setDofNr(dof); dof+=Ns;

  }
  us Piston::getNEqs() const {
    TRACE(15,"void Piston::getNEqs()");
    // 1 for equation of motion of piston
    // 2 for mass conservation left and right volumes
    // 2 for energy conservation left and right volumes
    // 2 for equations of state

    us neqs=7*Ns;
    // If not connected, one to set mass flow left and right to zero
    if(!leftConnected)
      neqs+=Ns;                 // For ml=0
    if(!rightConnected)
      neqs+=Ns;                 // For mr=0
    return neqs;
  }
  us Piston::getNDofs() const {
    TRACE(15,"us Piston::getNDofs()");
    return 10*Ns;
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

    rhor_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    Tr_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    pr_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;
    mr_.setadata(res.subvec(dofnr,dofnr+Ns-1)); dofnr+=Ns;

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

    res.subvec(dofnr,dofnr+Ns-1)=rhor_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=Tr_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=pr_(); dofnr+=Ns;
    res.subvec(dofnr,dofnr+Ns-1)=mr_(); dofnr+=Ns;

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
    const vd Vlt=V0l+xpt*Sl;
    const vd Vrt=V0r-xpt*Sr;

    // dVdt in time domain of left volume
    const vd dVldtt=iDFT*DDTfd*(xp_()*Sl);
    // dVdt in time domain of right volume
    const vd dVrdtt=iDFT*DDTfd*(-xp_()*Sr);

    // ************************************************************
    // Equation of motion of the piston
    error.subvec(eqnr,eqnr+Ns-1)=\
      M*DDTfd*DDTfd*xp_()\      // "Inertia" force
      +Cm*DDTfd*xp_()\          // Damping force
      +Km*xp_()\                // Spring force
      -pl_()*Sl\                // "Pushes"
      +pr_()*Sr\                // "Pushes back"
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
    }
    else{
      error.subvec(eqnr,eqnr+Ns-1)=DDTfd*fDFT*(rholt%Vlt)+ml_();
    }
    eqnr+=Ns;

    // ************************************************************
    // Mass conservation right side
    if(!rightConnected){
      assert(massR>0);
      error.subvec(eqnr,eqnr+Ns-1)=fDFT*(rhor_.tdata()%(V0r-xp_.tdata()*Sr));
      error(eqnr)-=massR;

      // Not connected right side, so set mr to zero
      eqnr+=Ns;
      error.subvec(eqnr,eqnr+Ns-1)=mr_();

    }
    else{
      error.subvec(eqnr,eqnr+Ns-1)=DDTfd*fDFT*(rhort%Vrt)+mr_();
    }
    eqnr+=Ns;

    // ************************************************************
    // Energy consrvation left side
    error.subvec(eqnr,eqnr+Ns-1)=gamfac*DDTfd*fDFT*(plt%Vlt)
      +fDFT*(plt%dVldtt)\ // Work contribution to change in
                                  \  // energy
      +fDFT*(cp*mlt%Tlt);            // Enthalpy flow out of segment

    if(!leftConnected){
      // Overwrite time-average part with constraint on time-average
      // internal energy, by fixing the time-averaged temperature
      error(eqnr)=Tl_()(0)-T0;

      // Isentropic model
      // error.subvec(eqnr,eqnr+Ns-1)=pl_()/p0;
      // error(eqnr)+=1;
      // error.subvec(eqnr,eqnr+Ns-1)+=-fDFT*pow(rhol_.tdata()/rho0,gamma);
    }
    eqnr+=Ns;

    // ************************************************************
    // Energy consrvation right side
    error.subvec(eqnr,eqnr+Ns-1)=gamfac*DDTfd*fDFT*(prt%Vrt)
      +fDFT*(prt%dVrdtt)\ // Work contribution to change in
                                  \  // energy
      +fDFT*(cp*mrt%Trt);            // Enthalpy flow out of segment

    if(!leftConnected){
      // Overwrite time-average part with constraint on time-average
      // internal energy, by fixing the time-averaged temperature
      error(eqnr)=Tr_()(0)-T0;

      // Isentropic model
      // error.subvec(eqnr,eqnr+Ns-1)=pr_()/p0;
      // error(eqnr)+=1;
      // error.subvec(eqnr,eqnr+Ns-1)+=-fDFT*pow(rhor_.tdata()/rho0,gamma);
    }
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
    const vd Vlt=V0l+xpt*Sl;
    const vd Vrt=V0r-xpt*Sr;

    // dVdt in time domain of left volume
    const vd dVldtt=iDFT*DDTfd*(xp_()*Sl);
    // dVdt in time domain of right volume
    const vd dVrdtt=iDFT*DDTfd*(-xp_()*Sr);

    // ************************************************************
    // 1: Equation of motion of the piston
    JacRow eomjac(eqnr,4);
    eomjac+=JacCol(xp_,M*DDTfd*DDTfd+Cm*DDTfd+Km*eye);
    eomjac+=JacCol(Fp_,-eye);
    eomjac+=JacCol(pl_,-Sl*eye);
    eomjac+=JacCol(pr_,Sr*eye);
    jac+=eomjac;
    eqnr+=Ns;

    // ************************************************************
    // 2: Mass conservation left side
    if(!leftConnected){
      JacRow mcljac(eqnr,2);
      mcljac+=JacCol(rhol_,fDFT*diagmat(V0l+xp_.tdata()*Sl)*iDFT);
      mcljac+=JacCol(xp_,fDFT*diagmat(rhol_.tdata()*Sl)*iDFT);
      jac+=mcljac;

      // Not connected to leftside so ml to zero
      eqnr+=Ns;
      jac+=JacRow(eqnr,JacCol(ml_,eye));
    }
    else{
      throw MyError("Not yet implemented");
    }
    eqnr+=Ns;

    // ************************************************************    
    // 3: Mass conservation right side
    if(!rightConnected){
      JacRow mcrjac(eqnr,2);
      mcrjac+=JacCol(rhor_,fDFT*diagmat(Vrt)*iDFT);
      mcrjac+=JacCol(xp_,fDFT*diagmat(-rhor_.tdata()*Sr)*iDFT);
      jac+=mcrjac;

      // Not connected to right side so mr is set to zero
      eqnr+=Ns;
      jac+=JacRow(eqnr,JacCol(mr_,eye));

    }
    else{
      throw MyError("Not yet implemented");
    }
    eqnr+=Ns;

    // ************************************************************
    // Energy equation left side:
    // 1/(gamma-1)*d/dt(p*Vl)+p*dVl/dt+ml*hl=0

    // 4: energy conservation left side
    // 4 columns in this row if not connected, otherwise 3
    JacRow enl(eqnr,2+(!leftConnected?1:0));
    
    // Terms due to change in internal energy
    dmat enlmat_pl=gamfac*DDTfd*fDFT*diagmat(Vlt)*iDFT;
    dmat enlmat_x=gamfac*DDTfd*fDFT*diagmat(plt*Sl)*iDFT;

    // Terms due to work done on fluid
    enlmat_pl+=fDFT*diagmat(dVldtt)*iDFT;
    enlmat_x+=fDFT*diagmat(Sl*plt)*iDFT*DDTfd;
    
    if(!leftConnected) {
      // Overwrite time-average part with constraint on time-average
      // internal energy, by fixing the time-averaged temperature
      dmat enlmat_T(Ns,Ns,fillwith::zeros);
      enlmat_T(0,0)=1;

      // Zero out first row
      enlmat_pl.row(0).zeros();
      enlmat_x.row(0).zeros();

      // Add Jacobian terms corresponding to left temperature
      enl+=JacCol(Tl_,enlmat_T);          

    }
    enl+=JacCol(pl_,enlmat_pl);
    enl+=JacCol(xp_,enlmat_x);    
    jac+=enl;
    eqnr+=Ns;

    // ************************************************************
    // Energy equation right side
    // 1/(gamma-1)*d/dt(p*Vr)+p*dVr/dt+mr*hr=0
    // 3 columns in this row if not connected, otherwise 2
    JacRow enr(eqnr,2+(!leftConnected?1:0));
    
    // Terms due to change in internal energy
    dmat enrmat_pr=gamfac*DDTfd*fDFT*diagmat(Vrt)*iDFT;
    dmat enrmat_x=gamfac*DDTfd*fDFT*diagmat(-prt*Sr)*iDFT;

    // Terms due to work done on fluid
    enrmat_pr+=fDFT*diagmat(dVrdtt)*iDFT;
    enrmat_x+=-fDFT*diagmat(Sr*prt)*iDFT*DDTfd;
    
    if(!rightConnected) {
      // Overwrite time-average part with constraint on time-average
      // internal energy, by fixing the time-averaged temperature
      dmat enrmat_T(Ns,Ns,fillwith::zeros);
      enrmat_T(0,0)=1;

      // Zero out first row
      enrmat_pr.row(0).zeros();
      enrmat_x.row(0).zeros();

      // Add Jacobian terms corresponding to left temperature
      enr+=JacCol(Tr_,enrmat_T);          

    }
    enr+=JacCol(pr_,enrmat_pr);
    enr+=JacCol(xp_,enrmat_x);    
    jac+=enr;
    eqnr+=Ns;

    const d Rs=gc->gas().Rs();
    // ************************************************************
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
  void Piston::show(us detailnr) const {
    TRACE(15,"void Piston::show()");
    cout << "Piston segment:\n"
         << "Piston mass: " << M
         << "\nPiston mechanical stiffness: "<< Km
         << "\nPiston mechanical damping: "<< Cm
         << "\nLeft volume: "<< V0l    
         << "\nRight volume: "<< V0r << endl;    
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
      WARN("Code has not been checked");
      mass+=arma::dot( fDFT.row(0) , rhol_.tdata()*(V0l+Sl*xp_.tdata()) );
    }
    if(rightConnected){
      WARN("Code has not been checked");
      mass+=arma::dot( fDFT.row(0) , rhor_.tdata()*(V0r-Sr*xp_.tdata()) );
    }
    return mass;
  }
  void Piston::dmtotdx(vd& dmtotdx) const {
    TRACE(15,"void Piston::dmtotdx()");
    // if(leftConnected){
    //   us rholdof=rhol.getDofNr();
    //   dmtotdx.subvec(rholdof,rholdof+Ns-1)=
    WARN("Not finished code");
    assert(false);
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
    return var(*gc,V0l+xp_.tdata()*Sl,false);
  }
  var Piston::Vr() const {
    TRACE(15,"var Piston::Vr()");
    return var(*gc,V0r-xp_.tdata()*Sr,false);
  }
} // namespace mech
//////////////////////////////////////////////////////////////////////
