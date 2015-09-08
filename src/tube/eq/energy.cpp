// energy.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

// #define TRACERPLUS 20
// #define ENERGY_SCALE (1/v.gc->rho0/v.gc->c0)
#define ENERGY_SCALE (1.0/v.gc->p0())
// #define ENERGY_SCALE (1.0/100)
// #define ENERGY_SCALE (1.0)
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

  vd kappaRt(const Cell& v) {		// Returns thermal conductivity time domain data
    TRACE(15,"Energy::kappaRt()");
    const WeightFactors w=WeightFactors(v);
    const vd& Tt=v.T().tdata();
    const vd& TtR=v.TR().tdata();
    return v.gc->gas().kappa(w.WrR*TtR+w.WrL*Tt);
  }
  vd kappaLt(const Cell& v) {		// Returns thermal conductivity time domain data
    TRACE(15,"Energy::kappaRt()");
    const WeightFactors w=WeightFactors(v);
    const vd& Tt=v.T().tdata();
    const vd& TtL=v.TL().tdata();
    return v.gc->gas().kappa(w.WlL*TtL+w.WlR*Tt);
  }

  void Energy::show() const{
    cout << "------------- Full energy equation\n";
    cout << "Wddt     :"<<Wddt<<"\n";
    // cout << "WkinLl   :"<<WkinLl<<"\n";
    // cout << "WkinLr   :"<<WkinLr<<"\n";
    // cout << "WkinRl   :"<< WkinRl<<"\n";
    // cout << "WkinRr   :"<< WkinRr<<"\n";
    // cout << "WclL     :"<<WclL<<"\n";
    // cout << "WclR     :"<<WclR<<"\n";
    // cout << "WcRl     :"<<WcRl<<"\n";
    // cout << "WcrR     :"<<WcrR<<"\n";

  }
  void Energy::init(){
    TRACE(8,"Energy::init(tube)");
    t=&v.getTube();

    d vx=v.vx;
    d xl=v.xl;
    d xr=v.xr;    

    Wddt=v.vVf;
    // VARTRACE(25,v.vVf);
    // Difference of factor half
    // ddtEkin=d/dt*(Vf*0.5*rho*u^2)=d/dt(1/Sf^2)*(Vf*0.5*rho*U^2)
    // =Sf*(1/Sf^2)*(Vf*0.5*mu)=d/dt*0.5*Vf/Sf*mu=d/dt*0.5*(xr-xl)*mu
    Wddtkin=0.5*(v.xr-v.xl);

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
    assert(v.gc);

    // Time derivative of static enthalpy in cell
    d gamma=this->gamma();
    // Time derivative of total energy
    vd error=(Wddt*DDTfd*v.p()())/(gamma-1.0)+\
      Wddtkin*DDTfd*v.mu()();
    // Total enthalpy flow in and out
    error+=mHr(v)-mHl(v);
    // Heat flow in and out
    error+=QR(v)-QL(v);
    // External heat    
    #ifndef NOHEAT
    assert(t);
    error-=(v.xr-v.xl)*t->heatSource().Qsf(v);
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
    domg_.subvec(dofnr,dofnr+Ns-1)=ENERGY_SCALE*
      ( Wddt/( (gamma-1.0)*v.gc->getomg() )*DDTfd*v.p()()+	\
		    Wddtkin*DDTfd*v.mu()()/v.gc->getomg() );                
  }
  JacRow Energy::jac() const{
    TRACE(6,"Energy::jac()");
    JacRow jac(dofnr,12);
    d gamma=this->gamma();
    // Time-derivative of internal energy
    jac+=JacCol(v.p(),(Wddt*DDTfd)/(gamma-1.0));
    jac+=JacCol(v.mu(),Wddtkin*DDTfd);

    // Enthalpy flow out minus in
    jac+=dmHr(v);
    jac+=(dmHl(v)*=-1);

    // Heat conduction part
    jac+=dQR(v);
    jac+=(dQL(v)*=-1);

    // Transverse heat transver
    #ifndef NOHEAT
    jac+=(t->heatSource().dQsf(v)*=-(v.xr-v.xl));
    #endif

    jac*=ENERGY_SCALE;
    return jac;    
  }
  
  inline d cp(const Cell& c) {
    return c.gc->gas().cp(c.gc->T0());
  }

  inline vd TmidLt(const Cell& v,const WeightFactors& w){
    return w.WlL*v.TL().tdata() + w.WlR*v.T().tdata();
  }
  inline vd TmidRt(const Cell& v,const WeightFactors& w){
    return w.WrR*v.TR().tdata() + w.WrL*v.T().tdata();
  }
  inline vd mEkinL(const Cell& v,const WeightFactors& w) {
    TRACE(15,"mEkinL()");
    const vd& rhoLt=v.rhoL().tdata();
    const vd& rhot=v.rho().tdata();
    d Sfl=v.Sfl;
    vd rholt=w.WlL*rhoLt+w.WlR*rhot;
    return fDFT*(pow(v.ml().tdata(),3)/pow(Sfl*rholt,2));
  }
  inline vd mEkinR(const Cell& v,const WeightFactors& w) {
    TRACE(15,"mEkinR()");
    const vd& rhoRt=v.right()->rho().tdata();
    const vd& rhot=v.rho().tdata();
    d Sfr=v.Sfr;
    vd rhort=w.WrR*rhoRt+w.WrL*rhot;
    return fDFT*(pow(v.mr().tdata(),3)/pow(Sfr*rhort,2));
  }  
  inline JacRow dmEkinL(const Cell& v,const WeightFactors& w){
    TRACE(15,"dmEkinL()");
    JacRow dmEkinL(-1,4);

    const vd& rhot=v.rho().tdata();
    const vd& rhoLt=v.rhoL().tdata();
    d Sfl=v.Sfl;
    JacRow dEkinL(-1,3);
    d Sflsq=pow(Sfl,2);
    vd rholt_sq=(pow((w.WlL*rhoLt+w.WlR*rhot),2));
    vd rholt_tr=pow(w.WlL*rhoLt+w.WlR*rhot,3);
    const vd& mlt=v.ml().tdata();
    dmEkinL+=JacCol(v.ml(),fDFT*diagmat(1.5*pow(mlt,2)/\
					(rholt_sq*Sflsq))*iDFT);
    
    dmEkinL+=JacCol(v.rho(),-fDFT*\
		    diagmat(w.WlR*pow(mlt,3)/	\
			    (rholt_tr*Sflsq))*iDFT);
    dmEkinL+=JacCol(v.rhoL(),-fDFT*\
		    diagmat(w.WlL*pow(mlt,3)/	\
			    (rholt_tr*Sflsq))*iDFT);
    
    return dmEkinL;
  }
  inline JacRow dmEkinR(const Cell& v,const WeightFactors& w){
    TRACE(15,"dmEkinR()");
    JacRow dmEkinR(-1,3);
    const vd& rhot=v.rho().tdata();
    const vd& rhoRt=v.rhoR().tdata();
    const vd& mrt=v.mr().tdata();
    d Sfr=v.Sfr;
    d Sfrsq=pow(Sfr,2);
    vd rhortsq=pow(w.WrR*rhoRt+w.WrL*rhot,2);
    vd rhort_tr=pow(w.WrR*rhoRt+w.WrL*rhot,3);      
    dmEkinR+=JacCol(v.mr(),fDFT*\
		    diagmat(1.5*pow(mrt,2)/(rhortsq*Sfrsq))*iDFT);

    dmEkinR+=JacCol(v.rho(),-fDFT*\
		    diagmat(w.WrL*pow(mrt,3)/(rhort_tr*Sfrsq))*iDFT);
    dmEkinR+=JacCol(v.rhoR(),-fDFT*\
		    diagmat(w.WrR*pow(mrt,3)/(rhort_tr*Sfrsq))*iDFT);
    return dmEkinR;
  }
  vd Energy::mHl(const Cell& v) {
    TRACE(15,"Energy::mHl()");

    if(v.left()) {
      const WeightFactors w=WeightFactors(v);
      // ******************** Static enthalpy part
      d cp0=cp(v);
      // TRACE(20,TmidLt(v,w));
      vd mHl=fDFT*(cp0*TmidLt(v,w)%v.ml().tdata());
      // ******************** Kinetic energy part
      mHl+=mEkinL(v,w);
      return mHl;
    }
    else
      return static_cast<const BcCell&>(v).mHbc()();
  }
  vd Energy::mHr(const Cell& v) {
    TRACE(15,"Energy::mHr()");
    if(v.right()){
      const WeightFactors w=WeightFactors(v);
      // ******************** Static enthalpy part    
      // TRACE(20,TmidRt(v,w));
      d cp0=cp(v);
      vd mHr=fDFT*(cp0*TmidRt(v,w)%v.mr().tdata());
      // ******************** Kinetic energy part
      mHr+=mEkinR(v,w);
    
      return mHr;
    }
    else
      return static_cast<const BcCell&>(v).mHbc()();
  }
  JacRow Energy::dmHl(const Cell& v){
    TRACE(15,"Energy::dmHl()");
    if(v.left()){
      JacRow dmHl(-1,3);
      const WeightFactors w=WeightFactors(v);
      // ******************** Static enthalpy part    
      // TRACE(20,TmidLt(v,w));      
      d cp0=cp(v);
      const vd& mltdata=v.ml().tdata();
      dmHl+=JacCol(v.ml(),fDFT*diagmat(cp0*TmidLt(v,w))*iDFT );
      dmHl+=JacCol(v.T(),fDFT*diagmat(cp0*w.WlR*mltdata)*iDFT );
      dmHl+=JacCol(v.TL(),fDFT*diagmat(cp0*w.WlL*mltdata)*iDFT );
      // ******************** Kinetic energy part
      dmHl+=dmEkinL(v,w);
      return dmHl;
    }
    else
      return JacRow(JacCol(static_cast<const BcCell&>(v).mHbc(),eye(v)));

  }
  JacRow Energy::dmHr(const Cell& v) {
    TRACE(15,"Energy::dmHr(const Cell& v)");

    if(v.right()){
      JacRow dmHr(-1,3);
      const WeightFactors w=WeightFactors(v);
      // ******************** Static enthalpy part    
      // TRACE(20,TmidRt(v,w));
      d cp0=cp(v);
      const vd& mrtdata=v.mr().tdata();
      dmHr+=JacCol(v.mr(),fDFT*diagmat(cp0*TmidRt(v,w))*iDFT);

      dmHr+=JacCol(v.T(),fDFT*diagmat(cp0*w.WrL*mrtdata)*iDFT);
      dmHr+=JacCol(v.TR(),fDFT*diagmat(cp0*w.WrR*mrtdata)*iDFT);
      // ******************** Kinetic energy part
      dmHr+=dmEkinR(v,w);
      return dmHr;
    }
    else
      return JacRow(JacCol(static_cast<const BcCell&>(v).mHbc(),eye(v)));
  }
  class heatW{
  public:
    d WclL=0,WclR=0,WcrL=0,WcrR=0; // Conduction weight factors    // d WkinLl=0,WkinLr=0,WkinRl=0,WkinRr=0;
    heatW(const Cell& v){
      d vx=v.vx;
      d xl=v.xl;
      d xr=v.xr;    
      if(v.left()){
        d vxim1=v.left()->vx;
        WclL=v.Sfl/(vx-vxim1);
      }
      else{
        WclL=v.Sfl/vx;
      }
      WclR=-WclL;
    
      if(v.right()){    
        d vxip1=v.right()->vx;
        WcrL=v.Sfr/(vxip1-vx);
      }
      else{
        WcrL=v.Sfr/(v.xr-vx);
      }
      WcrR=-WcrL;
    }
  };

  vd Energy::QL(const Cell& v) {
    TRACE(15,"Energy::QL()");
    const heatW w(v);
    const vd& Tt=v.T().tdata();
    const vd& Ttl=v.TL().tdata();
    // VARTRACE(15,fDFT*(kappaLt_()%(w.WclL*Ttl+w.WclR*Tt)));
    return fDFT*(kappaLt(v)%(w.WclL*Ttl+w.WclR*Tt));
  }
  JacRow Energy::dQL(const Cell& v) {
    TRACE(15,"Energy::dQL()");
    const heatW w(v);
    JacRow dQL(2);
    dQL+=JacCol(v.T(),fDFT*diagmat(w.WclR*kappaLt(v))*iDFT);
    dQL+=JacCol(v.TL(),fDFT*diagmat(w.WclL*kappaLt(v))*iDFT);
    return dQL;
  }
  vd Energy::QR(const Cell& v) {
    TRACE(15,"Energy::QR()");
    const heatW w(v);
    const vd& Tt=v.T().tdata();
    const vd& Ttr=v.TR().tdata();
    return fDFT*(kappaRt(v)%(w.WcrL*Tt+w.WcrR*Ttr));
  }
  JacRow Energy::dQR(const Cell& v) {
    TRACE(15,"Energy::dQR()");
    const heatW w(v);
    JacRow dQR(2);
    dQR+=JacCol(v.T(),fDFT*diagmat(w.WcrL*kappaRt(v))*iDFT);
    dQR+=JacCol(v.TR(),fDFT*diagmat(w.WcrR*kappaRt(v))*iDFT);
    return dQR;
  }
  vd Energy::extrapolateHeatFlow(const Cell& v) {
    TRACE(5,"Energy::extrapolateHeatFlow()");
    const heatW w(v);
    assert((!v.left() && v.right()) || (v.left() && !v.right()));
    if(!v.left()){
      return fDFT*(kappaLt(v)%(w.WclL*v.TL().tdata()+w.WclR*v.T().tdata()));
    }
    else {
      return fDFT*(kappaRt(v)%(w.WcrL*v.T().tdata()+w.WcrR*v.TR().tdata()));
    }
  }
  JacRow Energy::dExtrapolateHeatFlow(const Cell& v){
    TRACE(5,"Energy::dExtrapolateHeatFlow()");
    assert((!v.left() && v.right()) || (v.left() && !v.right()));    
    JacRow dQb(-1,2);
    const heatW w(v);
    // VARTRACE(30,kappaLt)      ;
    // VARTRACE(30,kappaRt);
    if(!v.left()){
      vd kappaLt_=kappaLt(v);
      dQb+=JacCol(v.T(),fDFT*diagmat(w.WclR*kappaLt_)*iDFT);
      dQb+=JacCol(v.TL(),fDFT*diagmat(w.WclL*kappaLt_)*iDFT);
      return dQb;
    }
    else {
      vd kappaRt_=kappaRt(v);
      dQb+=JacCol(v.T(),fDFT*diagmat(w.WcrL*kappaRt_)*iDFT);
      dQb+=JacCol(v.TR(),fDFT*diagmat(w.WcrR*kappaRt_)*iDFT);
      return dQb;
    }
  }
  d Energy::gamma() const {
    // d T0=v.T(0);
    d T0=v.gc->T0(); 
    return v.gc->gas().gamma(T0);
  }

  vd ExtrapolateEnthalpyFlow::extrapolateEnthalpyFlow(const BcCell& v) {
    TRACE(15,"ExtrapolateEnthalpyflow::extrapolateEnthalpyFlow()");
    assert((!v.left() && v.right()) || (!v.right() && v.left()));
    vd half_u_sqt=0.5*pow(v.ubc().tdata(),2);
    d cp0=cp(v);
    return fDFT*(v.mbc().tdata()%(cp0*v.Tbc().tdata()+half_u_sqt));
  }
  JacRow ExtrapolateEnthalpyFlow::dExtrapolateEnthalpyFlow(const BcCell& v) {
    TRACE(15,"ExtrapolateEnthalpyFlow::dExtrapolateEnthalpyFlow()");
    assert((!v.left() && v.right()) || (!v.right() && v.left()));
    JacRow jac(-1,4);
    const vd& mbct=v.mbc().tdata();
    const vd& ubct=v.ubc().tdata();
    const vd& Tbct=v.Tbc().tdata();
    d cp0=cp(v);
    d Sfbc=v.Sfbc();
    d Sfbcsq=pow(Sfbc,2);
    jac+=JacCol(v.mbc(),fDFT*diagmat(cp0*Tbct+0.5*pow(ubct,2))*iDFT);
    jac+=JacCol(v.Tbc(),fDFT*diagmat(cp0*mbct)*iDFT);
    jac+=JacCol(v.ubc(),fDFT*diagmat(mbct%ubct)*iDFT);
    return jac;
  }


} // namespace tube

  //////////////////////////////////////////////////////////////////////

