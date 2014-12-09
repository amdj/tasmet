// #define ENERGY_SCALE (1/v.gc->rho0/v.gc->c0)
// #define ENERGY_SCALE (1.0/v.gc->p0)
// #define ENERGY_SCALE (1.0/100)
#define ENERGY_SCALE (1.0)

#ifdef NOHEAT
#error Noheat already defined!
#endif
// #define NOHEAT

#include "energyeq.h"
#include "tubevertex.h"
#include "weightfactors.h"
#include "tube.h"
#include "jacobian.h"


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
    cout << "Wddt    :"<<Wddt<<"\n";
    cout << "Wgim1   :"<<Wgim1<<"\n";
    cout << "Wgim    :"<<Wgim<<"\n";
    cout << "Wgip    :"<<Wgip<<"\n";
    cout << "Wgip1   :"<<Wgip1<<"\n";
    cout << "Wkinim1 :"<<Wkinim1<<"\n";
    cout << "Wkini   :"<<Wkini<<"\n";
    cout << "Wkinip1 :"<<Wkinip1<<"\n";
    cout << "Wc1     :"<<Wc1<<"\n";
    cout << "Wc2     :"<<Wc2<<"\n";
    cout << "Wc3     :"<<Wc3<<"\n";
    cout << "Wc4     :"<<Wc4<<"\n";


  }
  void Energy::init(const WeightFactors& w,const Tube& t){
    TRACE(8,"Energy::init(tube)");
    heat=&t.getHeatSource();

    // d SfL=w.SfL;
    // d SfR=w.SfR;
    // d SfLsq=pow(SfL,2);
    // d SfRsq=pow(SfR,2);

    // const d& vSfR=rw.vSf;
    // const d& vSfL=lw.vSf;
    d vSfLsq=pow(w.vSfL,2);
    d vSfRsq=pow(w.vSfR,2);
    d vSfsq=pow(w.vSf,2);
    // d vSfLav=0.5*(lg.vSf+llg.vSf);
    // d vSfRav=0.5*(lg.vSf+rlg.vSf);


    Wddt=w.vVf;
    Wddtkin=0.5*Wddt/pow(w.vSf,2);

    Wgim1=-w.wLl;
    Wgim =-w.wLr;
    Wgip = w.wRl;
    Wgip1= w.wRr;

    Wkinim1=-0.5*w.wLl/vSfLsq;
    Wkini=0.5*(w.wRl/vSfsq-w.wLr/vSfsq);
    Wkinip1=0.5*w.wRr/vSfRsq;

    // Wkinim1=-0.5*w.wLl/SfLsq;
    // Wkini=0.5*(w.wRl/SfRsq-w.wLr/SfLsq);
    // Wkinip1=0.5*w.wRr/SfRsq;
    

    // TRACE(1,"dxm:"<< dxm);

    // TRACE(1,"dxp:"<< dxp);
    if(v.left()){
      WcLL=-w.SfL/(w.vx-w.vxm1);
      WcL=  w.SfL/(w.vx-w.vxm1);
    }
    else{
      WcLL=-w.SfL/(w.vx);
      WcL=  w.SfL/(w.vx);
    }
    if(left()&&right()){

      Wc2= w.SfL/(w.vx-w.vxm1);
      Wc3= w.SfR/(w.vxp1-w.vx);
      Wc4=-w.SfR/(w.vxp1-w.vx);
    }
    else if(v.left()){


}
  }
  JacRow Energy::jac() const{
    TRACE(6,"Energy::jac()");
    JacRow jac(dofnr,12);
    TRACE(0,"Energy, dofnr jac:"<< dofnr);
    jac+=drhoi();
    jac+=dTi();
    jac+=dpR();
    jac+=dpL();
    jac+=dUi();

    if(v.left()){
      jac+=drhoim1();
      jac+=dTim1();
      jac+=dUim1();
    }
    if(v.right()){
      jac+=drhoip1();
      jac+=dTip1();
      jac+=dUip1();
    }
    jac*=ENERGY_SCALE;
    return jac;    
  }
  vd Energy::QR() const{
    return kappaR()*(TR()()+T()())/

  }
  d Energy::Htot() const{
    TRACE(10,"Energy::Htot()");
    const dmat& fDFT=v.gc->fDFT;
    d Htot=0;
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    Htot+=gamfac*0.5*((v.pL()+v.pR())*v.U())(0);

    if(v.left()&&v.right()){
      const vd& Tim1=v.TL().tdata();
      const vd& Tip1=v.TR().tdata();
      vd conduction=fDFT*(0.5*v.localGeom().vSf*(kappaL()+kappaR())%(Tim1-Tip1))/(v.localGeom().xR+v.localGeom().xL);
      Htot+=conduction(0);
    }
    return Htot;
  }
  
  vd Energy::error() const {		// Error in momentum equation
    TRACE(6,"Energy::Error(), i="<<v.geti());
    assert(v.gc!=NULL);
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
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

    // Thermal energy flux
    error+=Wgim*fDFT*(gamfac*pLt%Uti);    
    error+=Wgip*fDFT*(gamfac*pRt%Uti);
    
    // Flux of kinetic energy
    error+=fDFT*(Wkini*rhoti%pow(Uti,3));

    // Conduction
    error+=fDFT*((Wc2*kappaL()+Wc3*kappaR())%Tti);

    if(v.left()){
      const vd& Utim1=v.UL().tdata();
      const vd& rhotim1=v.rhoL().tdata();

      // Thermal energy flux
      error+=fDFT*(Wgim1*gamfac*pLt%Utim1);
      // Specially for boundary conditions
      error+=fDFT*(WgUim1pR*gamfac*pRt%Utim1);
      
      // Flux of kinetic energy
      error+=fDFT*(Wkinim1*rhotim1%pow(Utim1,3));

      // Conduction term
      const vd& Tim1=v.TL().tdata();
      error+=Wc1*fDFT*(kappaL()%Tim1);

    }
    // TRACE(100,"error:"<<error);
    if(v.right()){
      const vd& Utip1=v.UR().tdata();
      const vd& rhotip1=v.rhoR().tdata(); 

      // Thermal energy flux
      error+=fDFT*(Wgip1*gamfac*pRt%Utip1);

      // Specially for boundary conditions
      error+=fDFT*(WgUip1pL*gamfac*pLt%Utip1);      
      
      // Flux of kinetic energy
      error+=fDFT*(Wkinip1*rhotip1%pow(Utip1,3));

      // Conduction term
      const vd& Tip1=v.TR().tdata();
      error+=Wc4*fDFT*(kappaR()%Tip1);
    }

    // Artificial viscosity terms      
    #ifdef EN_VISCOSITY
    TRACE(100,"Energy viscosity defined");
    if(v.left!=NULL && v.right!=NULL){
      error+=d_l()*(v.cWart2*v.p()+v.cWart1*v.left()->p() );
      error+=d_r()*(v.cWart3*v.p()+v.cWart4*v.right()->p() );
    }
    #endif

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
  void Energy::domg(vd& domg_) const {
    TRACE(0,"Energy::domg()");
    assert(v.gc!=NULL);
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
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
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      

    // Wddt = vVf, as defined in tubevertex.cpp
    d gamfac=gamma/(gamma-1.0);
    dpL+=(0.5*Wddt/(gamma-1.0))*DDTfd;
    dpL+=Wgim*gamfac*fDFT*v.U().diagt()*iDFT;
    if(v.left()){
      dpL+=(Wgim1*gamfac)*fDFT*v.UL().diagt()*iDFT;
    }

    // Specially for boundary conditions
    if(v.right())
      dpL+=(WgUip1pL*gamfac)*fDFT*v.UR().diagt()*iDFT;      

    return dpL;
  }
  JacCol Energy::dpR() const {
    TRACE(0,"Energy::dpR()");
    d gamma=this->gamma();
    JacCol dpR(v.pR());
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      

    d gamfac=gamma/(gamma-1.0);
    // Wddt = vVf, as defined in tubevertex.cpp

    dpR+=(0.5*Wddt/(gamma-1.0))*DDTfd;

    dpR+=Wgip*gamfac*fDFT*v.U().diagt()*iDFT;
    if(v.right()){
      dpR+=(Wgip1*gamfac)*fDFT*v.UR().diagt()*iDFT;
    }
    if(v.left()){
      dpR+=(WgUim1pR*gamfac)*fDFT*v.UL().diagt()*iDFT;
    }
    
    return dpR;
  }
  JacCol Energy::dUi() const {
    TRACE(0,"Energy::dUi()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
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
    dUi+=fDFT*diagmat(Wgim*gamfac*pLt)*iDFT;       
    dUi+=fDFT*diagmat(Wgip*gamfac*pRt)*iDFT; 

    // Flux of kinetic energy
    dUi+=fDFT*diagmat(Wkini*3*rhoti%pow(Uti,2))*iDFT;

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
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    JacCol drhoi(v.rho());			    // Initialize with zeros
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);

    const vd& Uti=v.U().tdata();
    drhoi+=DDTfd*fDFT*diagmat(Wddtkin*pow(Uti,2))*iDFT; // This term should get a eWuddt factor later on
    drhoi+=fDFT*diagmat(Wkini*pow(Uti,3))*iDFT;

    return drhoi;
  }  
  JacCol Energy::dUim1() const {
    TRACE(0,"Energy::dUim1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);

    const vd& rhotim1=v.rhoL().tdata(); 
    JacCol dUim1(v.UL());
    const vd& Utim1=v.UL().tdata();

    vd pLt=v.pL().tdata()+v.getp0t();
    vd pRt=v.pR().tdata()+v.getp0t();
    dUim1+=(Wgim1*gamfac)*fDFT*diagmat(pLt)*iDFT;

    // Specially for boundary conditions
    dUim1+=fDFT*diagmat(WgUim1pR*gamfac*pRt)*iDFT;    
    
    dUim1+=fDFT*diagmat(Wkinim1*3*rhotim1%pow(Utim1,2))*iDFT;

    
    return dUim1;
  }
  JacCol Energy::dUip1() const {
    TRACE(0,"Energy::dUip1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    JacCol dUip1(v.UR());
    const vd& Utip1=v.UR().tdata();
    vd pLt=v.pL().tdata()+v.getp0t();    
    vd pRt=v.pR().tdata()+v.getp0t();

    const vd& rhotip1=v.rhoR().tdata(); 

    dUip1+=fDFT*diagmat(Wgip1*gamfac*pRt)*iDFT;

    // Specially for boundary conditions
    dUip1+=fDFT*diagmat(WgUip1pL*gamfac*pLt)*iDFT;    
    
    dUip1+=fDFT*diagmat(Wkinip1*3*rhotip1%pow(Utip1,2))*iDFT;

    return dUip1;
  }
  JacCol Energy::drhoim1() const {
    TRACE(0,"Energy::drhoim1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    JacCol drhoim1(v.rhoL());
    const vd& Utim1=v.UL().tdata();
    drhoim1+=fDFT*diagmat(Wkinim1*pow(Utim1,3))*iDFT;
    return drhoim1;
  }
  JacCol Energy::drhoip1() const {
    TRACE(0,"Energy::drhoip1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    JacCol drhoip1(v.rhoR());

    const vd& rhotip1=v.rhoR().tdata();
    const vd& Utip1=v.UR().tdata();
    drhoip1+=fDFT*diagmat(Wkinip1*pow(Utip1,3))*iDFT;

    return drhoip1;
  }

  // ############################## TEMPERATURE AND CONDUCTION TERMS
  JacCol Energy::dTip1() const {
    TRACE(0,"Energy::dTip1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      

    JacCol dTip1(v.TR());
    // TRACE(100,"dTip1, right is"<<v.right);
    dTip1+=Wc4*fDFT*diagmat(kappaR())*iDFT;

    return dTip1;
  }
  JacCol Energy::dTi() const {
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    TRACE(0,"Energy::dTi()");
    JacCol dTi(v.T());
    dTi+=fDFT*diagmat(Wc2*kappaL()+Wc3*kappaR())*iDFT;
    assert(heat!=NULL);
    #ifndef NOHEAT
    dTi+=Wddt*heat->dTi(v);
    #endif
    return dTi;
  }
  JacCol Energy::dTim1() const {
    TRACE(0,"Energy::dTim1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    JacCol dTim1(v.TL());
    dTim1+=Wc1*fDFT*diagmat(kappaL())*iDFT;
    return dTim1;
  }

  // From here on, auxiliary functions
  vd Energy::kappaL()  const {
    TRACE(5,"Energy::kappaL()");
    // const dmat& DDTfd=v.gc->DDTfd;
    // const dmat& fDFT=v.gc->fDFT;
    // const dmat& iDFT=v.gc->iDFT;      
    
    vd kappaL(v.gc->Ns(),fillwith::zeros);
    const vd& Tti=v.T().tdata();
    vd kappait=v.gc->gas.kappa(Tti);

    if(v.left()){
      const vd& Ttim1=v.TL().tdata();
      vd kappaitm1=v.gc->gas.kappa(Ttim1);
      kappaL=wLr*kappait+wLl*kappaitm1;
    }
    else{
      const vd& Ttip1=v.TR().tdata();
      vd kappaitp1=v.gc->gas.kappa(Ttip1);
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
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    
    vd kappaR(v.gc->Ns(),fillwith::zeros);
    const vd& Tti=v.T().tdata();
    vd kappait=v.gc->gas.kappa(Tti);

    if(v.right()){
      const vd& Ttip1=v.TR().tdata();
      vd kappaitp1=v.gc->gas.kappa(Ttip1);
      kappaR=wRl*kappait+wRr*kappaitp1;
    }
    else{
      const vd& Ttim1=v.TL().tdata();
      vd kappaitm1=v.gc->gas.kappa(Ttim1);
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
