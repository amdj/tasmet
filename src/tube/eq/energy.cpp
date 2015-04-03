// #define TRACERPLUS 20
// #define ENERGY_SCALE (1/v.gc->rho0/v.gc->c0)
// #define ENERGY_SCALE (1.0/v.gc->p0)
// #define ENERGY_SCALE (1.0/100)
#define ENERGY_SCALE (1.0)

#ifdef NOHEAT
#error Noheat already defined!
#endif
#define NOHEAT

#include "cell.h"
#include "weightfactors.h"
#include "tube.h"
#include "jacrow.h"
#include "energy.h"

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)
#define DDTfd (v.gc->DDTfd)

namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  void Energy::show() const{
    cout << "------------- Full energy equation\n";
    cout << "Wddt     :"<<Wddt<<"\n";
    cout << "WLl     :"<<WLl<<"\n";
    cout << "WLr     :"<<WLr<<"\n";
    cout << "WRl     :"<<WRl<<"\n";
    cout << "WRr     :"<<WRr<<"\n";
    // cout << "WkinLl   :"<<WkinLl<<"\n";
    // cout << "WkinLr   :"<<WkinLr<<"\n";
    // cout << "WkinRl   :"<< WkinRl<<"\n";
    // cout << "WkinRr   :"<< WkinRr<<"\n";
    cout << "WcLl     :"<<WcLl<<"\n";
    cout << "WcLr     :"<<WcLr<<"\n";
    cout << "WcRl     :"<<WcRl<<"\n";
    cout << "WcRr     :"<<WcRr<<"\n";

  }
  void Energy::init(){
    TRACE(8,"Energy::init(tube)");
    const Tube& t=v.getTube();
    heat=&t.getHeatSource();

    d vx=v.vx;
    d xL=v.xL;
    d xR=v.xR;    

    Wddt=v.vVf;
    // VARTRACE(25,v.vVf);
    // Difference of factor half
    Wddtkin=0.5*(v.xR-v.xL);

    // WkinL=-0.5*w.wLl/SfLsq;
    // Wkin=0.5*(w.wRl/SfRsq-w.wLr/SfLsq);
    // WkinR=0.5*w.wRr/SfRsq;
   
    // vxm1=0 for the leftmost cell, so this is always true:
    // WcLl=w.SfL/(w.vx-w.vxm1);
    // WcLr=-WcLl;
    TRACE(8,"Energy::init(tube)");

    // WkinL=-0.5*w.wLl/SfLsq;
    // Wkin=0.5*(w.wRl/SfRsq-w.wLr/SfLsq);
    // WkinR=0.5*w.wRr/SfRsq;
    

    if(v.left()){
      d vxim1=v.left()->vx;
      WLl=(vx-xL)/(vx-vxim1);
      WLr=1-WLl;
      WcLl=v.SfL/(vx-vxim1);
      //   d vSfLsq=pow(w.vSfL,2);
      //   WkinLl=0.5*w.wLl/vSfLsq;
      //   WkinLr=0.5*w.wLr/vSfsq;
    }
    else{
      //   d SfLsq=pow(w.SfL,2);
      WcLl=v.SfL/vx;
      WcLr=-WcLl;
      //   WkinLl=0.5/SfLsq;
      //   WkinLr=0;
      WLl=1;
      WLr=0;
    }
    WcLr=-WcLl;
    
    if(v.right()){    
      d vxip1=v.right()->vx;
      WRr=(xR-vx)/(vxip1-vx);
      WRl=1-WRr;
      //   d vSfRsq=pow(w.vSfR,2);

      WcRl=v.SfR/(vxip1-vx);

      //   WkinRl=0.5*(w.wRl/vSfsq);
      //   WkinRr=0.5*(w.wRr/vSfRsq);

      //   WRl= w.wRl;
      //   WRr= w.wRr;
    }
    else{
    //   d SfRsq=pow(w.SfR,2);
      WcRl=v.SfR/(v.xR-vx);
    //   WkinRl=0;
    //   WkinRr=0.5/SfRsq;
      WRl=0;
      WRr=1;
    }
    WcRr=-WcRl;
  }
  d Energy::Htot() const{
    TRACE(10,"Energy::Htot()");
    // vd H=0.5*(HL()+HR());
    // return H(0);
    return 0;
  }
  vd Energy::error() const {		// Error in momentum equation
    TRACE(6,"Energy::Error(), i="<<v.geti());
    assert(v.gc!=nullptr);

    // Time derivative of static enthalpy in cell
    d gamma=this->gamma();
    vd error=(Wddt*DDTfd*v.p()())/(gamma-1.0);
    // Time derivative of total energy
    // error+=Wddtkin*DDTfd*v.mu()();

    // Total enthalpy flux
    error+=HR()-HL();
    // error+=hR()-hL();
    error+=QR()-QL();
    // VARTRACE(25,QR());
    // VARTRACE(25,QL());

    // External heat    
    assert(heat!=nullptr);
    #ifndef NOHEAT
    error+=Wddt*heat->heat(v);
    #else
    if(v.geti()==0)
      WARN("Applying no heat coupling");
    #endif

    // (Boundary source term)
    // error+=v.esource();
    return error;
  }
  JacRow Energy::jac() const{
    TRACE(6,"Energy::jac()");
    JacRow jac(dofnr,12);
    d gamma=this->gamma();
    // Time-derivative of internal energy
    jac+=JacCol(v.p(),(Wddt*DDTfd)/(gamma-1.0));
    jac+=JacCol(v.mu(),Wddtkin*DDTfd);

    // Enthalpy flow out minus in
    jac+=dHR();
    jac+=(dHL()*=-1);

    // Heat conduction part
    jac+=dQR();
    jac+=(dQL()*=-1);

    // Transverse heat transver
    #ifndef NOHEAT
    jac+=JacCol(v.U(),Wddt*heat->dUi(v));
    jac+=JacCol(v.T(),Wddt*heat->dTi(v));
    #endif

    // jac*=ENERGY_SCALE;
    return jac;    
  }
  namespace 
  {
    inline d cp(const Cell& c) {
      return c.gc->gas().cp(c.gc->T0());
    }
  } // namespace 
  
  vd Energy::HL() const{
    TRACE(2,"Energy::HL()");

    // ******************** Static enthalpy part
    d cp0=cp(v);
    // Time domain temperature at left cell wall
    vd dTmidt=WLl*v.TL().tdata() + WLr*v.T().tdata();
    vd HL=fDFT*(cp0*dTmidt%v.mL().tdata());
    // ******************** Kinetic energy part

    
    return HL;
  }
  vd Energy::HR() const{
    TRACE(2,"Energy::HR()");
    // ******************** Static enthalpy part    
    d cp0=cp(v);
    // Time domain temperature at left cell wall
    vd dTmidt=WRr*v.TR().tdata() + WRl*v.T().tdata();
    vd HR=fDFT*(cp0*dTmidt%v.mR().tdata());
    // ******************** Kinetic energy part

    // Time domain mu/(rho*Sf)
    // vd mu_ov_rhoSf_t=wRl*v.mu().tdata()/(v.SfR*v.rho().tdata())+
      // wRr*v.mu().tdata()/(v.SfR*v.rho().tdata())+
    
    return HR;
  }
  JacRow Energy::dHL() const{
    TRACE(5,"Energy::dHL()");

    JacRow dHL(8);
    d cp0=cp(v);
    vd dTmidt=WLl*v.TL().tdata()+
      WLr*v.T().tdata();

    // ******************** Static enthalpy part    
    const vd& mLtdata=v.mL().tdata();
    dHL+=JacCol(v.mL(),fDFT*diagmat(cp0*dTmidt)*iDFT );
    dHL+=JacCol(v.TL(),fDFT*diagmat(cp0*WLl*mLtdata)*iDFT );
    dHL+=JacCol(v.T(),fDFT*diagmat(cp0*WLr*mLtdata)*iDFT );
    // ******************** Kinetic energy part

    return dHL;
  }
  JacRow Energy::dHR() const{
    JacRow dHR(8);

    // ******************** Static enthalpy part    
    d cp0=cp(v);

    // Time domain temperature at left cell wall
    vd dTmidt=WRr*v.TR().tdata()
      +WRl*v.T().tdata();
    
    const vd& mRtdata=v.mR().tdata();
    dHR+=JacCol(v.mR(),fDFT*diagmat(cp0*dTmidt)*iDFT);
    dHR+=JacCol(v.TR(),fDFT*diagmat(cp0*WRr*mRtdata)*iDFT);
    dHR+=JacCol(v.T(),fDFT*diagmat(cp0*WRl*mRtdata)*iDFT);
    // ******************** Kinetic energy part

    return dHR;
  }
  vd Energy::QL() const{
    TRACE(4,"Energy::QL()");
    const vd& Tt=v.T().tdata();
    const vd& Ttl=v.TL().tdata();
    VARTRACE(15,fDFT*(kappaLt()%(WcLl*Ttl+WcLr*Tt)));
    return fDFT*(kappaLt()%(WcLl*Ttl+WcLr*Tt));
  }
  JacRow Energy::dQL() const{
    TRACE(4,"Energy::dQL()");
    JacRow dQL(2);
    dQL+=JacCol(v.T(),fDFT*diagmat(WcLr*kappaLt())*iDFT);
    dQL+=JacCol(v.TL(),fDFT*diagmat(WcLl*kappaLt())*iDFT);
    return dQL;
  }
  vd Energy::QR() const{
    TRACE(4,"Energy::QR()");
    const vd& Tt=v.T().tdata();
    const vd& Ttr=v.TR().tdata();
    VARTRACE(15,fDFT*(kappaRt()%(WcRl*Tt+WcRr*Ttr)));
    return fDFT*(kappaRt()%(WcRl*Tt+WcRr*Ttr));
  }
  JacRow Energy::dQR() const{
    TRACE(4,"Energy::dQR()");
    JacRow dQR(2);
    dQR+=JacCol(v.T(),fDFT*diagmat(WcRl*kappaRt())*iDFT);
    dQR+=JacCol(v.TR(),fDFT*diagmat(WcRr*kappaRt())*iDFT);
    return dQR;
  }
  vd Energy::kappaRt()  const {		// Returns thermal conductivity time domain data
    TRACE(5,"Energy::kappaRt()");
    const vd& Tt=v.T().tdata();
    const vd& TtR=v.TR().tdata();
    // VARTRACE(25,v.gc->gas().kappa(WRr*TtR+WRl*Tt));
    return v.gc->gas().kappa(WRr*TtR+WRl*Tt);
  }
  vd Energy::kappaLt()  const {		// Returns thermal conductivity time domain data
    TRACE(5,"Energy::kappaRt()");
    const vd& Tt=v.T().tdata();
    const vd& TtL=v.TL().tdata();    
    // VARTRACE(25,v.gc->gas().kappa(WLl*TtL+WLr*Tt));
    return v.gc->gas().kappa(WLl*TtL+WLr*Tt);
  }
  void Energy::domg(vd& domg_) const {
    TRACE(0,"Energy::domg()");
    assert(v.gc!=nullptr);

    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    WARN("Is this correct??");
    vd domg_full=Wddt*DDTfd*v.p()()/(gamma-1.0)/v.gc->getomg(); // Static
    domg_full+=Wddtkin*DDTfd*v.mu()()/v.gc->getomg();                                                        // enthalpy
                                                             // term
    // domg_.subvec(dofnr+1,dofnr+2)=domg_full.subvec(1,2);
    domg_.subvec(dofnr,dofnr+v.gc->Ns()-1)=domg_full;
  }
  vd Energy::extrapolateHeatFlow() const{
    TRACE(5,"Energy::extrapolateHeatFlow()");
    vd Qb(v.gc->Ns());
    if(!v.left()){
      vd kappaLt=this->kappaLt();
      Qb=fDFT*(kappaLt%(WcLl*v.TL().tdata()+WcLr*v.T().tdata()));
    }
    else if(!v.right()){
      vd kappaRt=this->kappaRt();
      Qb=fDFT*(kappaRt%(WcRl*v.T().tdata()+WcRr*v.TR().tdata()));
    }
    else{
      WARN("Extrapolation of heat flow asked for an inner cell. This is a bug. Aborting...");
      abort();
    }
    return Qb;
  }
  JacRow Energy::dExtrapolateHeatFlow() const{
    TRACE(5,"Energy::dExtrapolateHeatFlow()");
    JacRow dQb(2);

    // VARTRACE(30,kappaLt)      ;
    // VARTRACE(30,kappaRt);
    if(!v.left()){
      vd kappaLt=this->kappaLt();
      dQb+=JacCol(v.T(),fDFT*diagmat(WcLr*kappaLt)*iDFT);
      dQb+=JacCol(v.TL(),fDFT*diagmat(WcLl*kappaLt)*iDFT);
    }
    else if(!v.right()){
      vd kappaRt=this->kappaRt();
      dQb+=JacCol(v.T(),fDFT*diagmat(WcRl*kappaRt)*iDFT);
      dQb+=JacCol(v.TR(),fDFT*diagmat(WcRr*kappaRt)*iDFT);
    }
    else{
      WARN("That went fatally wrong!");
      abort();
    }
    return dQb;
  }
  d Energy::gamma() const {
    // d T0=v.T(0);
    d T0=v.gc->T0(); 
    return v.gc->gas().gamma(T0);
  }


} // namespace tube
