// #define ENERGY_SCALE (1/v.gc->rho0/v.gc->c0)
// #define ENERGY_SCALE (1.0/v.gc->p0)
// #define ENERGY_SCALE (1.0/100)
#define ENERGY_SCALE (1.0)

#ifdef NOHEAT
#error Noheat already defined!
#endif
// #define NOHEAT

#include "tubevertex.h"
#include "weightfactors.h"
#include "tube.h"
#include "jacobian.h"
#include "energy.h"

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)

namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  void Energy::show() const{
    cout << "------------- Full energy equation\n";
    #ifdef EN_VISCOSITY
    cout << "Artificial viscosity turned ON for energy equation\n";
    #else
    cout << "Artificial viscosity turned OFF for energy equation\n";
    #endif
    cout << "Wddt     :"<<Wddt<<"\n";
    cout << "WgLl     :"<<WgLl<<"\n";
    cout << "WgLr     :"<<WgLr<<"\n";
    cout << "WgRl     :"<<WgRl<<"\n";
    cout << "WgRr     :"<<WgRr<<"\n";
    cout << "WkinLl   :"<<WkinLl<<"\n";
    cout << "WkinLr   :"<<WkinLr<<"\n";
    cout << "WkinRl   :"<< WkinRl<<"\n";
    cout << "WkinRr   :"<< WkinRr<<"\n";
    cout << "WcLl     :"<<WcLl<<"\n";
    cout << "WcLr     :"<<WcLr<<"\n";
    cout << "WcRl     :"<<WcRl<<"\n";
    cout << "WcRr     :"<<WcRr<<"\n";


  }
  void Energy::init(){
    TRACE(8,"Energy::init(tube)");
    const Tube& t=v.getTube();
    heat=&t.getHeatSource();
    const WeightFactors& w=v.weightFactors();
    d vSfLsq=pow(w.vSfL,2);
    d vSfRsq=pow(w.vSfR,2);
    d vSfsq=pow(w.vSf,2);

    Wddt=w.vVf;
    Wddtkin=0.5*Wddt/pow(w.vSf,2);

    // WkinL=-0.5*w.wLl/SfLsq;
    // Wkin=0.5*(w.wRl/SfRsq-w.wLr/SfLsq);
    // WkinR=0.5*w.wRr/SfRsq;
    
    if(v.left()){
      WcLl=w.SfL/(w.vx-w.vxm1);
      WcLr=-w.SfL/(w.vx-w.vxm1);

      WkinLl=0.5*w.wLl/vSfLsq;
      WkinLr=0.5*w.wLr/vSfsq;

      WgLl=w.wLl;
      WgLr=w.wLr;

    }
    else{
      WcLl=w.SfL/(w.vx);
      WcLr=-w.SfL/(w.vx);

      WkinLl=-0.5/vSfLsq;
      WkinLr=0;

      WgLl=-1;
      WgLr=0;

    }
    if(v.right()){    
      WcRl= w.SfR/(w.vxp1-w.vx);
      WcRr=-w.SfR/(w.vxp1-w.vx);

      WkinRl=0.5*(w.wRl/vSfsq);
      WkinRr=0.5*w.wRr/vSfRsq;

      WgRl= w.wRl;
      WgRr= w.wRr;

    }
    else{
      WcRl= w.SfR/(w.xR-w.vx);
      WcRr=-w.SfR/(w.xR-w.vx);

      WkinRl=0;
      WkinRr=0.5/vSfRsq;

      WgRl=0;
      WgRr=1;
    }
  }
  JacRow Energy::jac() const{
    TRACE(6,"Energy::jac()");
    JacRow jac(dofnr,12);
    TRACE(0,"Energy, dofnr jac:"<< dofnr);

    jac+=(dHL()*=-1);
    jac+=dHR();
    
    jac+=drhoi();
    jac+=dTi();
    jac+=dpR();
    jac+=dpL();
    jac+=dUi();

    if(v.left()){
      jac+=drhoL();
      jac+=dTL();
      jac+=dUL();
    }
    if(v.right()){
      jac+=drhoR();
      jac+=dTR();
      jac+=dUR();
    }
    jac*=ENERGY_SCALE;
    return jac;    
  }
  // vd Energy::QR() const{
  //   return kappaR()%(TR()()+T()())/

  // }
  d Energy::Htot() const{
    TRACE(10,"Energy::Htot()");
    vd H=0.5*(HL()+HR());
    return H(0);
  }
  
  vd Energy::error() const {		// Error in momentum equation
    TRACE(6,"Energy::Error(), i="<<v.geti());
    assert(v.gc!=NULL);
    const dmat& DDTfd=v.gc->DDTfd;

    vd error(v.gc->Ns(),fillwith::zeros);
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);

    const vd& Uti=v.U().tdata();
    vd pti=0.5*(v.pL().tdata()+v.pR().tdata())+v.getp0t();
    vd pRt=v.getp0t()+v.pR().tdata();
    vd pLt=v.getp0t()+v.pL().tdata();
    const vd& Tti=v.T().tdata();
    const vd& rhoti=v.rho().tdata();

    // Time derivative of thermal energy
    error+=Wddt*DDTfd*0.5*(v.pL()()+v.pR()())/(gamma-1.0);

    // Time derivative of kinetic energy
    error+=Wddtkin*DDTfd*fDFT*(rhoti%Uti%Uti); //

    // Total enthalpy flux
    error+=HR()-HL();

    // External heat    
    assert(heat!=NULL);
    #ifndef NOHEAT
    error+=Wddt*heat->heat(v);
    #else
    if(v.i==0)
      TRACE(25,"Applying no heat coupling");
    #endif
    // (Boundary source term)
    error+=v.esource();
    // TRACE(100,"error:"<<error);
    return ENERGY_SCALE*error;
  }
  vd Energy::hL() const{
    TRACE(2,"Energy::hL()");
    // We still assume gamma is constant
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    const vd& UtL=v.UL().tdata();
    const vd& pLt=v.pL().tdata();
    const vd& Ut=v.U().tdata();
    return gamfac*fDFT*(pLt%(WgLl*UtL+WgLr*Ut));
    }
  JacRow Energy::dhL() const{
    JacRow dhL(3);
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    const vd& UtL=v.UL().tdata();
    const vd& pLt=v.pL().tdata();
    const vd& Ut=v.U().tdata();    
    dhL+=JacCol(v.pL(),fDFT*diagmat(WgLl*UtL+WgLr*Ut)*iDFT);
    dhL+=JacCol(v.U(),fDFT*(WgLr*v.pL().diagt())*iDFT);
    dhL+=JacCol(v.UL(),fDFT*(WgLl*v.pL().diagt())*iDFT);
    return dhL;
  }
  JacRow Energy::dHL() const{
    JacRow dHL=dhL();
    const vd& UL=v.UL().tdata();
    dHL+=JacCol(v.rhoL(),WkinLl*fDFT*(diagmat(UL%UL)*iDFT));
  }
  vd Energy::hR() const{
    TRACE(2,"Energy::hR()");
    // We still assume gamma is constant
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    const vd& UtR=v.UR().tdata();
    const vd& pRt=v.pR().tdata();
    const vd& Ut=v.rho().tdata();    
    return gamfac*fDFT*(pRt%(WgRl*Ut+WgRr*UtR));
  }
  vd Energy::HL() const{
    TRACE(2,"Energy::HL()");
    vd HL=hL();

    const vd& rhot=v.rho().tdata();
    const vd& Ut=v.rho().tdata();    
    const vd& rhotL=v.rhoL().tdata();
    const vd& UtL=v.UL().tdata();    

    HL+=fDFT*(WkinLl*rhotL%pow(UtL,3));
    HL+=fDFT*(WkinLr*rhot%pow(Ut,3));
    return HL;
  }
  vd Energy::HR() const{
    TRACE(2,"Energy::HR()");
    vd HR=hR();

    const vd& rhot=v.rho().tdata();
    const vd& Ut=v.rho().tdata();    
    const vd& rhotR=v.rhoR().tdata();
    const vd& UtR=v.UR().tdata();    

    HR+=fDFT*(WkinRl*rhot%pow(Ut,3));
    HR+=fDFT*(WkinRr*rhotR%pow(UtR,3));
    return HR;
  }
  vd Energy::EkinL() const{
    return HL()-hL();
  }
  vd Energy::EkinR() const{
    return HR()-hR();
  }

  void Energy::domg(vd& domg_) const {
    TRACE(0,"Energy::domg()");
    assert(v.gc!=NULL);
    const dmat& DDTfd=v.gc->DDTfd;

    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    WARN("Is this correct??");
    vd domg_full=0.5*Wddt*DDTfd*(v.pL()()+v.pR()())/(gamma-1.0)/v.gc->getomg(); // Static
    domg_full+=Wddtkin*DDTfd*(fDFT*(v.rho().tdata()%v.U().tdata()%v.U().tdata()))/v.gc->getomg();                                                        // enthalpy
                                                             // term
    // domg_.subvec(dofnr+1,dofnr+2)=domg_full.subvec(1,2);
    domg_.subvec(dofnr,dofnr+v.gc->Ns()-1)=domg_full;
  }

  JacCol Energy::dpL() const {
    TRACE(0,"Energy::dpL()");
    d gamma=this->gamma();
    JacCol dpL(v.pL());
    const dmat& DDTfd=v.gc->DDTfd;
    // Wddt = vVf, as defined in tubevertex.cpp
    d gamfac=gamma/(gamma-1.0);
    dpL+=(0.5*Wddt/(gamma-1.0))*DDTfd;
    dpL+=WgLr*gamfac*fDFT*v.U().diagt()*iDFT;

    dpL+=(WgLl*gamfac)*fDFT*v.UL().diagt()*iDFT;

    // Specially for boundary conditions
    dpL+=(WgURpL*gamfac)*fDFT*v.UR().diagt()*iDFT;      

    return dpL;
  }
  JacCol Energy::dpR() const {
    TRACE(0,"Energy::dpR()");
    d gamma=this->gamma();
    JacCol dpR(v.pR());
    const dmat& DDTfd=v.gc->DDTfd;

    d gamfac=gamma/(gamma-1.0);
    // Wddt = vVf, as defined in tubevertex.cpp

    dpR+=(0.5*Wddt/(gamma-1.0))*DDTfd;

    dpR+=WgRl*gamfac*fDFT*v.U().diagt()*iDFT;
    if(v.right()){
      dpR+=(WgRr*gamfac)*fDFT*v.UR().diagt()*iDFT;
    }
    if(v.left()){
      dpR+=(WgULpR*gamfac)*fDFT*v.UL().diagt()*iDFT;
    }
    
    return dpR;
  }
  JacCol Energy::dUi() const {
    TRACE(0,"Energy::dUi()");
    const dmat& DDTfd=v.gc->DDTfd;

    JacCol dUi(v.U());			    // Initialize with zeros
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);

    const vd& rhoti=v.rho().tdata();
    const vd& Uti=v.U().tdata();
    vd pRt=v.pR().tdata()+v.getp0t();
    vd pLt=v.pL().tdata()+v.getp0t();    


    // Derivative of kinetic energy
    dUi+=DDTfd*fDFT*diagmat(2*Wddtkin*rhoti%Uti)*iDFT;

    // Flux of thermal energy
    dUi+=fDFT*diagmat(WgLr*gamfac*pLt)*iDFT;       
    dUi+=fDFT*diagmat(WgRl*gamfac*pRt)*iDFT; 

    // Flux of kinetic energy
    dUi+=fDFT*diagmat(Wkin*3*rhoti%pow(Uti,2))*iDFT;

    // External heat
    assert(heat!=NULL);
    #ifdef NOHEAT
    dUi+=Wddt*heat->dUi(v);
    #endif
    return dUi;
  }
  JacCol Energy::drhoi() const {
    TRACE(0,"Energy::drhoi()");
    const dmat& DDTfd=v.gc->DDTfd;

    JacCol drhoi(v.rho());			    // Initialize with zeros
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);

    const vd& Uti=v.U().tdata();
    drhoi+=DDTfd*fDFT*diagmat(Wddtkin*pow(Uti,2))*iDFT; // This term should get a eWuddt factor later on
    drhoi+=fDFT*diagmat(Wkin*pow(Uti,3))*iDFT;

    return drhoi;
  }  
  JacCol Energy::dUL() const {
    TRACE(0,"Energy::dUL()");

    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);

    const vd& rhotL=v.rhoL().tdata(); 
    JacCol dUL(v.UL());
    const vd& UtL=v.UL().tdata();

    vd pLt=v.pL().tdata()+v.getp0t();
    vd pRt=v.pR().tdata()+v.getp0t();
    dUL+=(WgLl*gamfac)*fDFT*diagmat(pLt)*iDFT;

    // Specially for boundary conditions
    dUL+=fDFT*diagmat(WgULpR*gamfac*pRt)*iDFT;    
    
    dUL+=fDFT*diagmat(WkinL*3*rhotL%pow(UtL,2))*iDFT;

    
    return dUL;
  }
  JacCol Energy::dUR() const {
    TRACE(0,"Energy::dUR()");

    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    JacCol dUR(v.UR());
    const vd& UtR=v.UR().tdata();
    vd pLt=v.pL().tdata()+v.getp0t();    
    vd pRt=v.pR().tdata()+v.getp0t();

    const vd& rhotR=v.rhoR().tdata(); 

    dUR+=fDFT*diagmat(WgRr*gamfac*pRt)*iDFT;

    // Specially for boundary conditions
    dUR+=fDFT*diagmat(WgURpL*gamfac*pLt)*iDFT;    
    
    dUR+=fDFT*diagmat(WkinR*3*rhotR%pow(UtR,2))*iDFT;

    return dUR;
  }
  JacCol Energy::drhoL() const {
    TRACE(0,"Energy::drhoL()");

    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    JacCol drhoL(v.rhoL());
    const vd& UtL=v.UL().tdata();
    drhoL+=fDFT*diagmat(WkinL*pow(UtL,3))*iDFT;
    return drhoL;
  }
  JacCol Energy::drhoR() const {
    TRACE(0,"Energy::drhoR()");

    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    JacCol drhoR(v.rhoR());

    const vd& rhotR=v.rhoR().tdata();
    const vd& UtR=v.UR().tdata();
    drhoR+=fDFT*diagmat(WkinR*pow(UtR,3))*iDFT;

    return drhoR;
  }

  // ############################## TEMPERATURE AND CONDUCTION TERMS
  JacCol Energy::dTR() const {
    TRACE(0,"Energy::dTR()");
    const dmat& DDTfd=v.gc->DDTfd;

    JacCol dTR(v.TR());
    // TRACE(100,"dTR, right is"<<v.right);
    dTR+=WcRr*fDFT*diagmat(kappaR())*iDFT;

    return dTR;
  }
  JacCol Energy::dTi() const {
    const dmat& DDTfd=v.gc->DDTfd;

    TRACE(0,"Energy::dTi()");
    JacCol dTi(v.T());
    dTi+=fDFT*diagmat(WcLr*kappaL()+WcRl*kappaR())*iDFT;
    assert(heat!=NULL);
    #ifndef NOHEAT
    dTi+=Wddt*heat->dTi(v);
    #endif
    return dTi;
  }
  JacCol Energy::dTL() const {
    TRACE(0,"Energy::dTL()");
    const dmat& DDTfd=v.gc->DDTfd;

    JacCol dTL(v.TL());
    dTL+=WcLl*fDFT*diagmat(kappaL())*iDFT;
    return dTL;
  }

  // From here on, auxiliary functions
  vd Energy::kappaL()  const {
    TRACE(5,"Energy::kappaL()");
    // const dmat& DDTfd=v.gc->DDTfd;
    
    vd kappaL(v.gc->Ns(),fillwith::zeros);
    const vd& Tti=v.T().tdata();
    vd kappait=v.gc->gas.kappa(Tti);

    if(v.left()){
      const vd& TtL=v.TL().tdata();
      vd kappaitm1=v.gc->gas.kappa(TtL);
      kappaL=wLr*kappait+wLl*kappaitm1;
    }
    else{
      const vd& TtR=v.TR().tdata();
      vd kappaitp1=v.gc->gas.kappa(TtR);
      kappaL=wL1*kappaitp1+wL0*kappait;
    }
    // TRACE(100,"kappaL(0):"<<kappaL(0));
    // vd kappaL(v.gc->Ns(),fillwith::ones);
    // kappaL*=v.gc->gas.kappa(v.gc->T0);
    // kappaL.zeros();		// WARNING
    return kappaL;
  }
  vd Energy::kappaR()  const {		// Returns thermal conductivity time domain data
    TRACE(5,"Energy::kappaR()");
    const dmat& DDTfd=v.gc->DDTfd;
    
    vd kappaR(v.gc->Ns(),fillwith::zeros);
    const vd& Tti=v.T().tdata();
    vd kappait=v.gc->gas.kappa(Tti);

    if(v.right()){
      const vd& TtR=v.TR().tdata();
      vd kappaitp1=v.gc->gas.kappa(TtR);
      kappaR=wRl*kappait+wRr*kappaitp1;
    }
    else{
      const vd& TtL=v.TL().tdata();
      vd kappaitm1=v.gc->gas.kappa(TtL);
      kappaR=wRNm2*kappaitm1+wRNm1*kappait;
    }
    // TRACE(100,"kappaR(0):"<<kappaR(0));
    // vd kappaR(v.gc->Ns(),fillwith::ones);
    // kappaR*=v.gc->gas.kappa(v.gc->T0);

    // kappaR.zeros();		// WARNING!!
    return kappaR;
  }
  d Energy::gamma() const {
    // d T0=v.T(0);
    d T0=v.gc->T0; 
    return v.gc->gas.gamma(T0);
  }


} // namespace tube
