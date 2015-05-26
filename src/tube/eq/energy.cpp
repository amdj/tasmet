// energy.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

// #define TRACERPLUS 20
// #define ENERGY_SCALE (1/v.gc->rho0/v.gc->c0)
// #define ENERGY_SCALE (1.0/v.gc->p0())
// #define ENERGY_SCALE (1.0/100)
#define ENERGY_SCALE (1.0)
// #define TRACERPLUS 15
#ifdef NOHEAT
#warn Heat turned off!!
#endif

#include "energy.h"
#include "bccell.h"
#include "weightfactors.h"
#include "tube.h"
#include "jacrow.h"


#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)
#define DDTfd (v.gc->DDTfd)
#define Ns (v.gc->Ns())

namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  using std::ignore;

  void Energy::show() const{
    cout << "------------- Full energy equation\n";
    cout << "Wddt     :"<<Wddt<<"\n";
    // cout << "WLl     :"<<WLl<<"\n";
    // cout << "WLr     :"<<WLr<<"\n";
    // cout << "WRl     :"<<WRl<<"\n";
    // cout << "WRr     :"<<WRr<<"\n";
    // cout << "WkinLl   :"<<WkinLl<<"\n";
    // cout << "WkinLr   :"<<WkinLr<<"\n";
    // cout << "WkinRl   :"<< WkinRl<<"\n";
    // cout << "WkinRr   :"<< WkinRr<<"\n";
    // cout << "WcLl     :"<<WcLl<<"\n";
    // cout << "WcLr     :"<<WcLr<<"\n";
    // cout << "WcRl     :"<<WcRl<<"\n";
    // cout << "WcRr     :"<<WcRr<<"\n";

  }
  void Energy::init(){
    TRACE(8,"Energy::init(tube)");
    t=&v.getTube();

    d vx=v.vx;
    d xL=v.xL;
    d xR=v.xR;    

    Wddt=v.vVf;
    WheatQ=(xR-xL);
    // VARTRACE(25,v.vVf);
    // Difference of factor half
    Wddtkin=0.5*(v.xR-v.xL);

    TRACE(8,"Energy::init(tube)");
  }
  // d Energy::Htot() const{
  //   TRACE(10,"Energy::Htot()");
  //   // vd H=0.5*(HL()+HR());
  //   // return H(0);
  //   return 0;
  // }
  vd Energy::error() const {		// Error in momentum equation
    TRACE(15,"Energy::Error(), i="<<v.geti());
    assert(v.gc!=nullptr);

    // Time derivative of static enthalpy in cell
    d gamma=this->gamma();
    vd error=(Wddt*DDTfd*v.p()())/(gamma-1.0);
    // Time derivative of total energy
    error+=Wddtkin*DDTfd*v.mu()();
    // Total enthalpy flow in and out
    error+=mHR(v)-mHL(v);
    // Heat flow in and out
    error+=QR(v)-QL(v);
    // External heat    
    #ifndef NOHEAT
    assert(t);
    error-=t->getHeatSource().heat(v);
    #else
    if(v.geti()==0)
      WARN("Applying no heat coupling");
    #endif
    // (Boundary source term)
    // error+=v.esource();
    return error*ENERGY_SCALE;
  }
  void Energy::domg(vd& domg_) const {
    TRACE(18,"Energy::domg()");
    assert(v.gc!=nullptr);
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    vd domg_full=ENERGY_SCALE*(Wddt/((gamma-1.0)*v.gc->getomg()))*DDTfd*v.p()();
    domg_full+=Wddtkin*DDTfd*v.mu()()/v.gc->getomg();                
    domg_.subvec(dofnr,dofnr+Ns-1)=domg_full;
  }
  JacRow Energy::jac() const{
    TRACE(6,"Energy::jac()");
    JacRow jac(dofnr,12);
    d gamma=this->gamma();
    // Time-derivative of internal energy
    jac+=JacCol(v.p(),(Wddt*DDTfd)/(gamma-1.0));
    jac+=JacCol(v.mu(),Wddtkin*DDTfd);

    // Enthalpy flow out minus in
    jac+=dmHR(v);
    jac+=(dmHL(v)*=-1);

    // Heat conduction part
    jac+=dQR(v);
    jac+=(dQL(v)*=-1);

    // Transverse heat transver
    #ifndef NOHEAT
    jac+=JacCol(v.mL(),-0.5*Wddt*t->getHeatSource().dmi(v));
    jac+=JacCol(v.mR(),-0.5*Wddt*t->getHeatSource().dmi(v));
    jac+=JacCol(v.T(),-t->getHeatSource().dTi(v));
    #endif

    jac*=ENERGY_SCALE;
    return jac;    
  }
  
  inline d cp(const Cell& c) {
    return c.gc->gas().cp(c.gc->T0());
  }

  inline vd TmidLt(const Cell& v,const WeightFactors& w){
    return w.WLl*v.TL().tdata() + w.WLr*v.T().tdata();
  }
  inline vd TmidRt(const Cell& v,const WeightFactors& w){
    return w.WRr*v.TR().tdata() + w.WRl*v.T().tdata();
  }
  inline vd EkinLt(const Cell& v,const WeightFactors& w) {
    TRACE(15,"EkinLt(v,w)");
    assert(v.left());
    const vd& rhot=v.rho().tdata();
    d vSf=v.vSf;
    const vd& rhoLt=v.rhoL().tdata();
    d vSfL=v.left()->vSf;
    return 0.5*pow(
                   v.mL().tdata()/(w.WLl*(rhoLt*vSfL)+w.WLr*(rhot*vSf))
                   ,2);
  }
  inline vd EkinRt(const Cell& v,const WeightFactors& w) {
    TRACE(15,"EkinRt(v,w)");
    // returns 0.5*u^2 right
    assert(v.right());
    const vd& rhoRt=v.right()->rho().tdata();
    const vd& rhot=v.rho().tdata();
    d vSfR=v.right()->vSf;
    d vSf=v.vSf;
    return 0.5*pow(
                   v.mR().tdata()/(w.WRr*(rhoRt*vSfR)+w.WRl*(rhot*vSf))
                   ,2);
  }
  inline vd mEkinL(const Cell& v,const WeightFactors& w) {
    TRACE(15,"mEkinL()");
    return fDFT*(v.mL().tdata()%EkinLt(v,w));
  }
  inline vd mEkinR(const Cell& v,const WeightFactors& w) {
    TRACE(15,"mEkinR()");
    return fDFT*(v.mR().tdata()%EkinRt(v,w));
  }  
  inline JacRow dmEkinL(const Cell& v,const WeightFactors& w){
    TRACE(15,"dmEkinL()");
    JacRow dmEkinL(-1,4);

    const vd& rhot=v.rho().tdata();
    d vSf=v.vSf;
    const vd& rhoLt=v.rhoL().tdata();
    d vSfL=v.left()->vSf;
    JacRow dEkinL(-1,1);
    vd one_ov_rhoSfLt_sq=1/pow(w.WLl*(rhoLt*vSfL)+w.WLr*(rhot*vSf),2);
    vd one_ov_rhoSfLt_tr=1/pow(w.WLl*(rhoLt*vSfL)+w.WLr*(rhot*vSf),3);

    dmEkinL+=JacCol(v.mL(),fDFT*diagmat(1.5*pow(v.mL().tdata(),2)%one_ov_rhoSfLt_sq)*iDFT);
    // dmEkinL+=JacCol(v.rho(),-fDFT*diagmat(vSf*w.WLr*pow(v.mL().tdata(),3)%one_ov_rhoSfLt_tr)*iDFT);
    // dmEkinL+=JacCol(v.rhoL(),-fDFT*diagmat(vSfL*w.WLl*pow(v.mL().tdata(),3)%one_ov_rhoSfLt_tr)*iDFT);
    
    return dmEkinL;
  }
  inline JacRow dmEkinR(const Cell& v,const WeightFactors& w){
    TRACE(15,"dmEkinR()");
    JacRow dmEkinR(-1,3);
    const vd& rhot=v.rho().tdata();
    d vSf=v.vSf;
    const vd& rhoRt=v.rhoR().tdata();
    d vSfR=v.right()->vSf;
    vd one_ov_rhoSfRt_sq=1/pow(w.WRr*(rhoRt*vSfR)+w.WRl*(rhot*vSf),2);
    vd one_ov_rhoSfRt_tr=1/pow(w.WRr*(rhoRt*vSfR)+w.WRl*(rhot*vSf),3);      
    dmEkinR+=JacCol(v.mR(),fDFT*diagmat(1.5*pow(v.mR().tdata(),2)%one_ov_rhoSfRt_sq)*iDFT);

    // dmEkinR+=JacCol(v.rho(),-fDFT*diagmat(vSf*w.WRl*pow(v.mR().tdata(),3)%one_ov_rhoSfRt_tr)*iDFT);
    // dmEkinR+=JacCol(v.rhoR(),-fDFT*diagmat(vSfR*w.WRr*pow(v.mR().tdata(),3)%one_ov_rhoSfRt_tr)*iDFT);
      

    return dmEkinR;
  }
  vd Energy::mHL(const Cell& v) {
    TRACE(15,"Energy::mHL()");

    if(v.left()) {
      const WeightFactors w=WeightFactors(v);
      // ******************** Static enthalpy part
      d cp0=cp(v);
      // TRACE(20,TmidLt(v,w));
      vd mHL=fDFT*(cp0*TmidLt(v,w)%v.mL().tdata());
      // ******************** Kinetic energy part
      // mHL+=mEkinL(v,w);
      return mHL;
    }
    else
      return static_cast<const BcCell&>(v).mHbc()();
  }
  vd Energy::mHR(const Cell& v) {
    TRACE(15,"Energy::mHR()");
    if(v.right()){
      const WeightFactors w=WeightFactors(v);
      // ******************** Static enthalpy part    
      // TRACE(20,TmidRt(v,w));
      d cp0=cp(v);
      vd mHR=fDFT*(cp0*TmidRt(v,w)%v.mR().tdata());
      // ******************** Kinetic energy part
      // mHR+=mEkinR(v,w);
    
      return mHR;
    }
    else
      return static_cast<const BcCell&>(v).mHbc()();
  }
  JacRow Energy::dmHL(const Cell& v){
    TRACE(15,"Energy::dmHL()");
    if(v.left()){
      JacRow dmHL(-1,3);
      const WeightFactors w=WeightFactors(v);
      // ******************** Static enthalpy part    
      // TRACE(20,TmidLt(v,w));      
      d cp0=cp(v);
      const vd& mLtdata=v.mL().tdata();
      dmHL+=JacCol(v.mL(),fDFT*diagmat(cp0*TmidLt(v,w))*iDFT );
      dmHL+=JacCol(v.T(),fDFT*diagmat(cp0*w.WLr*mLtdata)*iDFT );
      dmHL+=JacCol(v.TL(),fDFT*diagmat(cp0*w.WLl*mLtdata)*iDFT );
      // ******************** Kinetic energy part
      // dmHL+=dmEkinL(v,w);
      return dmHL;
    }
    else
      return JacRow(JacCol(static_cast<const BcCell&>(v).mHbc(),eye(v)));

  }
  JacRow Energy::dmHR(const Cell& v) {
    TRACE(15,"Energy::dmHR(const Cell& v)");

    if(v.right()){
      JacRow dmHR(-1,3);
      const WeightFactors w=WeightFactors(v);
      // ******************** Static enthalpy part    
      // TRACE(20,TmidRt(v,w));
      d cp0=cp(v);
      const vd& mRtdata=v.mR().tdata();
      dmHR+=JacCol(v.mR(),fDFT*diagmat(cp0*TmidRt(v,w))*iDFT);

      dmHR+=JacCol(v.T(),fDFT*diagmat(cp0*w.WRl*mRtdata)*iDFT);
      dmHR+=JacCol(v.TR(),fDFT*diagmat(cp0*w.WRr*mRtdata)*iDFT);
      // ******************** Kinetic energy part
      // dmHR+=dmEkinR(v,w);
      return dmHR;
    }
    else
      return JacRow(JacCol(static_cast<const BcCell&>(v).mHbc(),eye(v)));
  }
  class heatW{
  public:
    d WcLl=0,WcLr=0,WcRl=0,WcRr=0; // Conduction weight factors    // d WkinLl=0,WkinLr=0,WkinRl=0,WkinRr=0;
    heatW(const Cell& v){
      d vx=v.vx;
      d xL=v.xL;
      d xR=v.xR;    
      if(v.left()){
        d vxim1=v.left()->vx;
        WcLl=v.SfL/(vx-vxim1);
      }
      else{
        WcLl=v.SfL/vx;
      }
      WcLr=-WcLl;
    
      if(v.right()){    
        d vxip1=v.right()->vx;
        WcRl=v.SfR/(vxip1-vx);
      }
      else{
        WcRl=v.SfR/(v.xR-vx);
      }
      WcRr=-WcRl;
    }
  };

  vd Energy::QL(const Cell& v) {
    TRACE(15,"Energy::QL()");
    const heatW w(v);
    const vd& Tt=v.T().tdata();
    const vd& Ttl=v.TL().tdata();
    // VARTRACE(15,fDFT*(kappaLt()%(w.WcLl*Ttl+w.WcLr*Tt)));
    return fDFT*(kappaLt(v)%(w.WcLl*Ttl+w.WcLr*Tt));
  }
  JacRow Energy::dQL(const Cell& v) {
    TRACE(15,"Energy::dQL()");
    const heatW w(v);
    JacRow dQL(2);
    dQL+=JacCol(v.T(),fDFT*diagmat(w.WcLr*kappaLt(v))*iDFT);
    dQL+=JacCol(v.TL(),fDFT*diagmat(w.WcLl*kappaLt(v))*iDFT);
    return dQL;
  }
  vd Energy::QR(const Cell& v) {
    TRACE(15,"Energy::QR()");
    const heatW w(v);
    const vd& Tt=v.T().tdata();
    const vd& Ttr=v.TR().tdata();
    return fDFT*(kappaRt(v)%(w.WcRl*Tt+w.WcRr*Ttr));
  }
  JacRow Energy::dQR(const Cell& v) {
    TRACE(15,"Energy::dQR()");
    const heatW w(v);
    JacRow dQR(2);
    dQR+=JacCol(v.T(),fDFT*diagmat(w.WcRl*kappaRt(v))*iDFT);
    dQR+=JacCol(v.TR(),fDFT*diagmat(w.WcRr*kappaRt(v))*iDFT);
    return dQR;
  }
  vd Energy::kappaRt(const Cell& v) {		// Returns thermal conductivity time domain data
    TRACE(15,"Energy::kappaRt()");
    const WeightFactors w=WeightFactors(v);
    const vd& Tt=v.T().tdata();
    const vd& TtR=v.TR().tdata();
    return v.gc->gas().kappa(w.WRr*TtR+w.WRl*Tt);
    // return ones(Ns)*v.gc->gas().kappa(v.gc->T0());
  }
  vd Energy::kappaLt(const Cell& v) {		// Returns thermal conductivity time domain data
    TRACE(15,"Energy::kappaRt()");
    const WeightFactors w=WeightFactors(v);
    const vd& Tt=v.T().tdata();
    const vd& TtL=v.TL().tdata();    
    return v.gc->gas().kappa(w.WLl*TtL+w.WLr*Tt);
    // return ones(Ns)*v.gc->gas().kappa(v.gc->T0());
  }
  vd Energy::extrapolateEnthalpyFlow(const Cell& v) {
    TRACE(15,"Energy::extrapolateEnthalpyFlow()");
    // Can only be called for leftmost or rightmost node
    assert((!v.left() && v.right()) || (v.left() && !v.right()));
    if(!v.left()){
      d W0,W1; std::tie(W0,W1)=BcWeightFactorsW(v);
      return W0*mHR(v)+W1*mHR(*v.right());
    }
    else{
      d WR1,WR2; std::tie(WR1,WR2)=BcWeightFactorsW(v);
      return WR1*mHL(v)+WR2*mHL(*v.left());
    }
  }
  JacRow Energy::dExtrapolateEnthalpyFlow(const Cell& v) {
    TRACE(15,"Energy::dExtrapolateEnthalpyFlow()");
    JacRow jac(2);
    if(!v.left()){
      d W0,W1; std::tie(W0,W1)=BcWeightFactorsW(v);
      jac+=(dmHR(v)*=W0);
      jac+=(dmHR(*v.right())*=W1);
    }
    else{
      d WR1,WR2; std::tie(WR1,WR2)=BcWeightFactorsW(v);
      jac+=(dmHL(v)*=WR1);
      jac+=(dmHL(*v.left())*=WR2);
    }
    return jac;
  }
  vd Energy::extrapolateHeatFlow(const Cell& v) {
    TRACE(5,"Energy::extrapolateHeatFlow()");
    const heatW w(v);
    vd Qb(Ns);
    if(!v.left()){
      vd kappaLt=Energy::kappaLt(v);
      Qb=fDFT*(kappaLt%(w.WcLl*v.TL().tdata()+w.WcLr*v.T().tdata()));
    }
    else if(!v.right()){
      vd kappaRt=Energy::kappaRt(v);
      Qb=fDFT*(kappaRt%(w.WcRl*v.T().tdata()+w.WcRr*v.TR().tdata()));
    }
    else{
      WARN("Extrapolation of heat flow asked for an inner cell. This is a bug. Aborting...");
      abort();
    }
    return Qb;
  }
  JacRow Energy::dExtrapolateHeatFlow(const Cell& v){
    TRACE(5,"Energy::dExtrapolateHeatFlow()");
    JacRow dQb(2);
    const heatW w(v);
    // VARTRACE(30,kappaLt)      ;
    // VARTRACE(30,kappaRt);
    if(!v.left()){
      vd kappaLt=Energy::kappaLt(v);
      dQb+=JacCol(v.T(),fDFT*diagmat(w.WcLr*kappaLt)*iDFT);
      dQb+=JacCol(v.TL(),fDFT*diagmat(w.WcLl*kappaLt)*iDFT);
    }
    else if(!v.right()){
      vd kappaRt=Energy::kappaRt(v);
      dQb+=JacCol(v.T(),fDFT*diagmat(w.WcRl*kappaRt)*iDFT);
      dQb+=JacCol(v.TR(),fDFT*diagmat(w.WcRr*kappaRt)*iDFT);
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

  //////////////////////////////////////////////////////////////////////

