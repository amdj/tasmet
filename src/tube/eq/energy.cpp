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
    cout << "WLl     :"<<WLl<<"\n";
    cout << "WLr     :"<<WLr<<"\n";
    cout << "WRl     :"<<WRl<<"\n";
    cout << "WRr     :"<<WRr<<"\n";
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

      WLl=w.wLl;
      WLr=w.wLr;

    }
    else{
      WcLl=w.SfL/(w.vx);
      WcLr=-w.SfL/(w.vx);

      WkinLl=-0.5/vSfLsq;
      WkinLr=0;

      WLl=-1;
      WLr=0;

    }
    if(v.right()){    
      WcRl= w.SfR/(w.vxp1-w.vx);
      WcRr=-w.SfR/(w.vxp1-w.vx);

      WkinRl=0.5*(w.wRl/vSfsq);
      WkinRr=0.5*w.wRr/vSfRsq;

      WRl= w.wRl;
      WRr= w.wRr;

    }
    else{
      WcRl= w.SfR/(w.xR-w.vx);
      WcRr=-w.SfR/(w.xR-w.vx);

      WkinRl=0;
      WkinRr=0.5/vSfRsq;

      WRl=0;
      WRr=1;
    }
  }
  d Energy::Htot() const{
    TRACE(10,"Energy::Htot()");
    vd H=0.5*(HL()+HR());
    return H(0);
  }
  
  vd Energy::error() const {		// Error in momentum equation
    TRACE(6,"Energy::Error(), i="<<v.geti());
    assert(v.gc!=NULL);

    vd error=ddtEtot();         // Time derivative of total energy
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

    return ENERGY_SCALE*error;
  }
  JacRow Energy::jac() const{
    TRACE(6,"Energy::jac()");
    JacRow jac(dofnr,12);
    TRACE(0,"Energy, dofnr jac:"<< dofnr);

    jac+=dddtEtot();
    jac+=(dHL()*=-1);
    jac+=dHR();

    jac+=JacCol(v.U(),Wddt*heat->dUi(v));
    jac+=JacCol(v.T(),Wddt*heat->dTi(v));
    jac*=ENERGY_SCALE;
    return jac;    
  }
  vd Energy::ddtEtherm() const {
    TRACE(2,"Energy::ddtEstat()");
    const dmat& DDTfd=v.gc->DDTfd;
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    return Wddt*DDTfd*0.5*(v.pL()()+v.pR()())/(gamma-1.0);
  }
  JacRow Energy::dddtEtherm() const {
    TRACE(2,"Energy::dddtEtherm()");
    const dmat& DDTfd=v.gc->DDTfd;
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    JacRow dddtEtherm(4);
    dddtEtherm+=JacCol(v.pL(),0.5*Wddt*DDTfd/(gamma-1.0));
    dddtEtherm+=JacCol(v.pR(),0.5*Wddt*DDTfd/(gamma-1.0));
    return dddtEtherm;
  }
  vd Energy::ddtEtot() const{
    TRACE(2,"Energy::ddtEtot()");
    const dmat& DDTfd=v.gc->DDTfd;
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);

    const vd& rhot=v.rho().tdata();
    const vd& Ut=v.U().tdata();
    // Time derivative of kinetic energy
    vd ddtEtot=ddtEtherm();
    ddtEtot+=Wddtkin*DDTfd*fDFT*(rhot%Ut%Ut); //
    return ddtEtot;
  }
  JacRow Energy::dddtEtot() const {
    TRACE(2,"Energy::dddtEtot()");
    const dmat& DDTfd=v.gc->DDTfd;
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    const vd& rhot=v.rho().tdata();
    const vd& Ut=v.U().tdata();
    JacRow dddtEtot=dddtEtherm();
    dddtEtot+=JacCol(v.U(),2.0*Wddtkin*DDTfd*fDFT*diagmat(rhot%Ut)*iDFT);
    dddtEtot+=JacCol(v.rho(),Wddtkin*DDTfd*fDFT*diagmat(Ut%Ut)*iDFT);  
    return dddtEtot;  
  }
  vd Energy::hL() const{
    TRACE(2,"Energy::hL()");
    // We still assume gamma is constant
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    const vd& UtL=v.UL().tdata();
    const vd& pLt=v.pL().tdata();
    const vd& Ut=v.U().tdata();
    return gamfac*fDFT*(pLt%(WLl*UtL+WLr*Ut));
    }
  JacRow Energy::dhL() const{
    JacRow dhL(3);
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    const vd& UtL=v.UL().tdata();
    const vd& pLt=v.pL().tdata();
    const vd& Ut=v.U().tdata();    
    dhL+=JacCol(v.pL(),fDFT*diagmat(WLl*UtL+WLr*Ut)*iDFT);
    dhL+=JacCol(v.U(),fDFT*(WLr*v.pL().diagt())*iDFT);
    dhL+=JacCol(v.UL(),fDFT*(WLl*v.pL().diagt())*iDFT);
    return dhL;
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
  JacRow Energy::dHL() const{
    JacRow dHL=dhL();
    const vd& rhot=v.rho().tdata();
    const vd& Ut=v.rho().tdata();    
    const vd& rhotL=v.rhoL().tdata();
    const vd& UtL=v.UL().tdata();    
    dHL+=JacCol(v.U(),2.0*WkinLr*fDFT*(diagmat(rhot%Ut%Ut)*iDFT));
    dHL+=JacCol(v.rho(),WkinLr*fDFT*(diagmat(Ut%Ut%Ut)*iDFT));

    dHL+=JacCol(v.rhoL(),WkinLl*fDFT*(diagmat(pow(UtL,3))*iDFT));
    dHL+=JacCol(v.UL(),3.0*WkinLl*fDFT*(diagmat(rhotL%UtL%UtL)*iDFT));
    return dHL;
  }
  vd Energy::hR() const{
    TRACE(2,"Energy::hR()");
    // We still assume gamma is constant
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    const vd& UtR=v.UR().tdata();
    const vd& pRt=v.pR().tdata();
    const vd& Ut=v.rho().tdata();    
    return gamfac*fDFT*(pRt%(WRl*Ut+WRr*UtR));
  }
  JacRow Energy::dhR() const{
    JacRow dhR(3);
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    const vd& UtR=v.UR().tdata();
    const vd& pRt=v.pR().tdata();
    const vd& Ut=v.U().tdata();    
    dhR+=JacCol(v.pR(),fDFT*diagmat(WRr*UtR+WRl*Ut)*iDFT);
    dhR+=JacCol(v.U(),fDFT*(WRl*v.pR().diagt())*iDFT);
    dhR+=JacCol(v.UR(),fDFT*(WRr*v.pR().diagt())*iDFT);
    return dhR;
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
  JacRow Energy::dHR() const{
    JacRow dHR=dhR();
    const vd& rhot=v.rho().tdata();
    const vd& Ut=v.rho().tdata();    
    const vd& rhotR=v.rhoR().tdata();
    const vd& UtR=v.UR().tdata();    
    dHR+=JacCol(v.U(),2.0*WkinRl*fDFT*(diagmat(rhot%Ut%Ut)*iDFT));
    dHR+=JacCol(v.rho(),WkinRl*fDFT*(diagmat(Ut%Ut%Ut)*iDFT));

    dHR+=JacCol(v.rhoR(),WkinRr*fDFT*(diagmat(pow(UtR,3))*iDFT));
    dHR+=JacCol(v.UR(),3.0*WkinRr*fDFT*(diagmat(rhotR%UtR%UtR)*iDFT));
    return dHR;
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

  // JacCol Energy::dpL() const {
  //   TRACE(0,"Energy::dpL()");
  //   d gamma=this->gamma();
  //   JacCol dpL(v.pL());
  //   const dmat& DDTfd=v.gc->DDTfd;
  //   // Wddt = vVf, as defined in tubevertex.cpp
  //   d gamfac=gamma/(gamma-1.0);
  //   dpL+=(0.5*Wddt/(gamma-1.0))*DDTfd;
  //   dpL+=WLr*gamfac*fDFT*v.U().diagt()*iDFT;

  //   // dpL+=(WLl*gamfac)*fDFT*v.UL().diagt()*iDFT;

  //   // Specially for boundary conditions
  //   // dpL+=(WURpL*gamfac)*fDFT*v.UR().diagt()*iDFT;      

  //   return dpL;
  // }
  // JacCol Energy::dpR() const {
  //   TRACE(0,"Energy::dpR()");
  //   d gamma=this->gamma();
  //   JacCol dpR(v.pR());
  //   const dmat& DDTfd=v.gc->DDTfd;

  //   d gamfac=gamma/(gamma-1.0);
  //   // Wddt = vVf, as defined in tubevertex.cpp

  //   dpR+=(0.5*Wddt/(gamma-1.0))*DDTfd;

  //   dpR+=WRl*gamfac*fDFT*v.U().diagt()*iDFT;
  //   if(v.right()){
  //     dpR+=(WRr*gamfac)*fDFT*v.UR().diagt()*iDFT;
  //   }
  //   if(v.left()){
  //     dpR+=(WULpR*gamfac)*fDFT*v.UL().diagt()*iDFT;
  //   }
    
  //   return dpR;
  // }
  // JacCol Energy::dUi() const {
  //   TRACE(0,"Energy::dUi()");
  //   const dmat& DDTfd=v.gc->DDTfd;

  //   JacCol dUi(v.U());			    // Initialize with zeros
  //   d gamma=this->gamma();
  //   d gamfac=gamma/(gamma-1.0);

  //   const vd& rhoti=v.rho().tdata();
  //   const vd& Uti=v.U().tdata();
  //   vd pRt=v.pR().tdata()+v.getp0t();
  //   vd pLt=v.pL().tdata()+v.getp0t();    


  //   // Derivative of kinetic energy
  //   dUi+=DDTfd*fDFT*diagmat(2*Wddtkin*rhoti%Uti)*iDFT;

  //   // Flux of thermal energy
  //   dUi+=fDFT*diagmat(WLr*gamfac*pLt)*iDFT;       
  //   dUi+=fDFT*diagmat(WRl*gamfac*pRt)*iDFT; 

  //   // Flux of kinetic energy
  //   dUi+=fDFT*diagmat(Wkin*3*rhoti%pow(Uti,2))*iDFT;

  //   // External heat
  //   assert(heat!=NULL);
  //   #ifdef NOHEAT
  //   dUi+=Wddt*heat->dUi(v);
  //   #endif
  //   return dUi;
  // }
  // JacCol Energy::drhoi() const {
  //   TRACE(0,"Energy::drhoi()");
  //   const dmat& DDTfd=v.gc->DDTfd;

  //   JacCol drhoi(v.rho());			    // Initialize with zeros
  //   d gamma=this->gamma();
  //   d gamfac=gamma/(gamma-1.0);

  //   const vd& Uti=v.U().tdata();
  //   drhoi+=DDTfd*fDFT*diagmat(Wddtkin*pow(Uti,2))*iDFT; // This term should get a eWuddt factor later on
  //   drhoi+=fDFT*diagmat(Wkin*pow(Uti,3))*iDFT;

  //   return drhoi;
  // }  
  // JacCol Energy::dUL() const {
  //   TRACE(0,"Energy::dUL()");

  //   d gamma=this->gamma();
  //   d gamfac=gamma/(gamma-1.0);

  //   const vd& rhotL=v.rhoL().tdata(); 
  //   JacCol dUL(v.UL());
  //   const vd& UtL=v.UL().tdata();

  //   vd pLt=v.pL().tdata()+v.getp0t();
  //   vd pRt=v.pR().tdata()+v.getp0t();
  //   dUL+=(WLl*gamfac)*fDFT*diagmat(pLt)*iDFT;

  //   // Specially for boundary conditions
  //   dUL+=fDFT*diagmat(WULpR*gamfac*pRt)*iDFT;    
    
  //   dUL+=fDFT*diagmat(WkinL*3*rhotL%pow(UtL,2))*iDFT;

    
  //   return dUL;
  // }
  // JacCol Energy::dUR() const {
  //   TRACE(0,"Energy::dUR()");

  //   d gamma=this->gamma();
  //   d gamfac=gamma/(gamma-1.0);
  //   JacCol dUR(v.UR());
  //   const vd& UtR=v.UR().tdata();
  //   vd pLt=v.pL().tdata()+v.getp0t();    
  //   vd pRt=v.pR().tdata()+v.getp0t();

  //   const vd& rhotR=v.rhoR().tdata(); 

  //   dUR+=fDFT*diagmat(WRr*gamfac*pRt)*iDFT;

  //   // Specially for boundary conditions
  //   dUR+=fDFT*diagmat(WURpL*gamfac*pLt)*iDFT;    
    
  //   dUR+=fDFT*diagmat(WkinR*3*rhotR%pow(UtR,2))*iDFT;

  //   return dUR;
  // }

  vd Energy::kappaR()  const {		// Returns thermal conductivity time domain data
    TRACE(5,"Energy::kappaR()");
    const dmat& DDTfd=v.gc->DDTfd;
    const vd& Tt=v.T().tdata();
    const vd& TtR=v.TR().tdata();    
    return v.gc->gas.kappa(WRr*TtR+WRl*Tt);
  }
  vd Energy::kappaL()  const {		// Returns thermal conductivity time domain data
    TRACE(5,"Energy::kappaR()");
    const dmat& DDTfd=v.gc->DDTfd;
    const vd& Tt=v.T().tdata();
    const vd& TtL=v.TL().tdata();    
    return v.gc->gas.kappa(WLl*TtL+WLr*Tt);
  }
  d Energy::gamma() const {
    // d T0=v.T(0);
    d T0=v.gc->T0; 
    return v.gc->gas.gamma(T0);
  }


} // namespace tube
