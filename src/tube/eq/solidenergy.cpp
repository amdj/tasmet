#include "solidenergy.h"
#include "jacrow.h"
#include "cell.h"
#include "tube.h"
#include "solid.h"
#include "weightfactors.h"

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)
#define DDTfd (v.gc->DDTfd)
#define Ns (v.gc->Ns())


namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  SolidEnergy::SolidEnergy(const Cell& v,const solids::Solid* s,d ksfrac):
      Equation(v),
      ksfrac(ksfrac),
      solid(s)
    {
      TRACE(15,"SolidEnergy::SolidEnergy()");
      if(v.vSs<=0){
	WARN("In Tube " << v.getTube().getName() << ", with ID " << v.getTube().getID() << ":\n");
	WARN("Solid has no cross-sectional area. Geometry problem?");
	throw MyError("Solid has no cross-sectional area. Geometry problem?");
      }
    }

  void SolidEnergy::init() {
    t=&v.getTube();
  }
  vd SolidEnergy::kappaRt(const Cell& v) const {		// Returns thermal conductivity time domain data
    TRACE(15,"Energy::kappaRt()");
    assert(solid);
    const WeightFactors w(v);
    const vd& Tt=v.Ts().tdata();
    const vd& TtR=v.TsR().tdata();
    VARTRACE(10,solid->kappa(w.WrR*TtR+w.WrL*Tt));
    return solid->kappa(w.WrR*TtR+w.WrL*Tt);
  }
  vd SolidEnergy::kappaLt(const Cell& v) const {		// Returns thermal conductivity time domain data
    TRACE(15,"Energy::kappaRt()");
    const WeightFactors w(v);
    const vd& Tt=v.Ts().tdata();
    const vd& TtL=v.TsL().tdata();
    VARTRACE(10,solid->kappa(w.WlL*TtL+w.WlR*Tt));
    return solid->kappa(w.WlL*TtL+w.WlR*Tt);
  }
  struct heatW{
    d WclL=0,WclR=0,WcrL=0,WcrR=0; // Conduction weight factors
    heatW(const Cell& v){
      d vx=v.vx;
      d xl=v.xl;
      d xr=v.xr;    
      if(v.left()){
        d vxim1=v.left()->vx;
        WclL=v.Ssl/(vx-vxim1);
      }
      else{
        WclL=v.Ssl/vx;
      }
      WclR=-WclL;
    
      if(v.right()){    
        d vxip1=v.right()->vx;
        WcrL=v.Ssr/(vxip1-vx);
      }
      else{
        WcrL=v.Ssr/(v.xr-vx);
      }
      WcrR=-WcrL;
    }
  };
  void SolidEnergy::show() const{
    cout << "------------- SolidEnergy equation\n";
    // cout << "WclL     :"<<WclL<<"\n";
    // cout << "WclR     :"<<WclR<<"\n";
    // cout << "WcRl     :"<<WcRl<<"\n";
    // cout << "WcrR     :"<<WcrR<<"\n";

  }

  vd SolidEnergy::error() const {		// Error in momentum equation
    TRACE(15,"SolidEnergy::Error(), i="<<v.geti());
    d rhoc=solid->rho(v.Ts()(0))*solid->c(v.Ts()(0));
    VARTRACE(10,v.vVs*rhoc*DDTfd*v.Ts()()+QR(v)-QL(v));
    return v.vVs*rhoc*DDTfd*v.Ts()()+QR(v)-QL(v) \
      +(v.xr-v.xl)*t->heatSource().Qsf(v);
  }
  JacRow SolidEnergy::jac() const {
    TRACE(15,"SolidEnergy::jac()");
    JacRow jac(dofnr,5);
    d rhoc=solid->rho(v.Ts()(0))*solid->c(v.Ts()(0));
    VARTRACE(40,v.vVs);
    jac+=JacCol(v.Ts(),v.vVs*rhoc*DDTfd);
    jac+=dQR(v);
    jac+=(dQL(v)*=-1);
    jac+=(t->heatSource().dQsf(v)*=(v.xr-v.xl));

    return jac;
  }
  
  void SolidEnergy::domg(vd& domg) const {
    TRACE(15,"SolidEnergy::domg()");
    d rhoc=solid->rho(v.Ts()(0))*solid->c(v.Ts()(0));
    domg.subvec(dofnr,dofnr+Ns-1)=(v.vVs*rhoc/v.gc->getomg())*DDTfd*v.Ts()();
  }
  vd SolidEnergy::QR(const Cell& v) const {
    TRACE(15,"SolidEnergy::QR()");
    const heatW w(v);
    const vd& Tt=v.Ts().tdata();
    const vd& Ttr=v.TsR().tdata();
    return fDFT*(ksfrac*kappaRt(v)%(w.WcrL*Tt+w.WcrR*Ttr));
  }
  vd SolidEnergy::QL(const Cell& v) const {
    TRACE(15,"SolidEnergy::QL()");
    const heatW w(v);
    const vd& Tt=v.Ts().tdata();
    const vd& Ttl=v.TsL().tdata();
    return fDFT*(kappaLt(v)%(w.WclL*Ttl+w.WclR*Tt));
  }
  JacRow SolidEnergy::dQL(const Cell& v) const {
    TRACE(15,"SolidEnergy::dQL()");
    const heatW w(v);
    JacRow dQL(-1,2);
    dQL+=JacCol(v.Ts(),fDFT*diagmat(w.WclR*kappaLt(v))*iDFT);
    dQL+=JacCol(v.TsL(),fDFT*diagmat(w.WclL*kappaLt(v))*iDFT);
    return dQL;
  }
  JacRow SolidEnergy::dQR(const Cell& v) const {
    TRACE(15,"SolidEnergy::dQR()");
    const heatW w(v);
    JacRow dQR(-1,2);
    dQR+=JacCol(v.Ts(),fDFT*diagmat(w.WcrL*kappaRt(v))*iDFT);
    dQR+=JacCol(v.TsR(),fDFT*diagmat(w.WcrR*kappaRt(v))*iDFT);
    return dQR;
  }
} // namespace tube


