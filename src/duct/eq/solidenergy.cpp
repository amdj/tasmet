#include "solidenergy.h"
#include "jacrow.h"
#include "cell.h"
#include "duct.h"
#include "solid.h"
#include "weightfactors.h"

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)
#define DDTfd (v.gc->DDTfd)
#define Ns (v.gc->Ns())
#define ENERGY_SCALE (1.0/(solid->rho(v.gc->T0())*solid->c(v.gc->T0())))

namespace duct{
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
	WARN("In Duct " << v.getDuct().getName() << ", with ID " << v.getDuct().getID() << ":\n");
	WARN("Solid has no cross-sectional area. Geometry problem?");
	throw MyError("Solid has no cross-sectional area. Geometry problem?");
      }
    }

  void SolidEnergy::init() {
    t=&v.getDuct();
  }
  vd SolidEnergy::kappaRt() const {		// Returns thermal conductivity time domain data
    TRACE(15,"Energy::kappaRt()");
    assert(solid);
    const WeightFactors w(v);
    const vd& Tt=v.Ts().tdata();
    const vd& TtR=v.TsR().tdata();
    VARTRACE(10,solid->kappa(w.WrR*TtR+w.WrL*Tt));
    return solid->kappa(w.WrR*TtR+w.WrL*Tt);
  }
  vd SolidEnergy::kappaLt() const {		// Returns thermal conductivity time domain data
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
    VARTRACE(10,v.vVs*rhoc*DDTfd*v.Ts()()+QR()-QL());
    vd vdeeye(Ns,fillwith::zeros); vdeeye(0)=1;
    return (v.vVs*rhoc*DDTfd*v.Ts()()+QR()-QL()		\
	    +(v.xr-v.xl)*t->heatSource().Qsf(v)-Qin*vdeeye)*ENERGY_SCALE;
  }
  JacRow SolidEnergy::jac() const {
    TRACE(15,"SolidEnergy::jac()");
    JacRow jac(dofnr,5);
    d rhoc=solid->rho(v.Ts()(0))*solid->c(v.Ts()(0));
    VARTRACE(40,v.vVs);
    jac+=JacCol(v.Ts(),v.vVs*rhoc*DDTfd);
    jac+=dQR();
    jac+=(dQL()*=-1);
    jac+=(t->heatSource().dQsf(v)*=(v.xr-v.xl));
    jac*=ENERGY_SCALE;
    return jac;
  }
  
  void SolidEnergy::domg(vd& domg) const {
    TRACE(15,"SolidEnergy::domg()");
    d rhoc=solid->rho(v.Ts()(0))*solid->c(v.Ts()(0));
    domg.subvec(dofnr,dofnr+Ns-1)=ENERGY_SCALE*\
      ((v.vVs*rhoc/v.gc->getomg())*DDTfd*v.Ts()());
  }
  vd SolidEnergy::QR() const {
    TRACE(15,"SolidEnergy::QR()");
    const heatW w(v);
    const vd& Tt=v.Ts().tdata();
    const vd& Ttr=v.TsR().tdata();
    return fDFT*(ksfrac*kappaRt()%(w.WcrL*Tt+w.WcrR*Ttr));
  }
  vd SolidEnergy::QL() const {
    TRACE(15,"SolidEnergy::QL()");
    const heatW w(v);
    const vd& Tt=v.Ts().tdata();
    const vd& Ttl=v.TsL().tdata();
    return fDFT*(kappaLt()%(w.WclL*Ttl+w.WclR*Tt));
  }
  JacRow SolidEnergy::dQL() const {
    TRACE(15,"SolidEnergy::dQL()");
    const heatW w(v);
    JacRow dQL(-1,2);
    dQL+=JacCol(v.Ts(),fDFT*diagmat(w.WclR*kappaLt())*iDFT);
    dQL+=JacCol(v.TsL(),fDFT*diagmat(w.WclL*kappaLt())*iDFT);
    return dQL;
  }
  JacRow SolidEnergy::dQR() const {
    TRACE(15,"SolidEnergy::dQR()");
    const heatW w(v);
    JacRow dQR(-1,2);
    dQR+=JacCol(v.Ts(),fDFT*diagmat(w.WcrL*kappaRt())*iDFT);
    dQR+=JacCol(v.TsR(),fDFT*diagmat(w.WcrR*kappaRt())*iDFT);
    return dQR;
  }
  vd SolidEnergy::extrapolateHeatFlow() const {
    TRACE(5,"Energy::extrapolateHeatFlow()");
    const heatW w(v);
    assert((!v.left() && v.right()) || (v.left() && !v.right()));
    if(!v.left()){
      return fDFT*(ksfrac*kappaLt()%(w.WclL*v.TsL().tdata()+w.WclR*v.Ts().tdata()));
    }
    else {
      return fDFT*(ksfrac*kappaRt()%(w.WcrL*v.Ts().tdata()+w.WcrR*v.TsR().tdata()));
    }
  }
  JacRow SolidEnergy::dExtrapolateHeatFlow() const {
    TRACE(5,"Energy::dExtrapolateHeatFlow()");
    assert((!v.left() && v.right()) || (v.left() && !v.right()));    
    JacRow dQb(-1,2);
    const heatW w(v);
    if(!v.left()){
      vd kappaLt_=kappaLt();
      dQb+=JacCol(v.Ts(),fDFT*diagmat(ksfrac*w.WclR*kappaLt_)*iDFT);
      dQb+=JacCol(v.TsL(),fDFT*diagmat(ksfrac*w.WclL*kappaLt_)*iDFT);
      return dQb;
    }
    else {
      vd kappaRt_=kappaRt();
      dQb+=JacCol(v.Ts(),fDFT*diagmat(ksfrac*w.WcrL*kappaRt_)*iDFT);
      dQb+=JacCol(v.TsR(),fDFT*diagmat(ksfrac*w.WcrR*kappaRt_)*iDFT);
      return dQb;
    }
  }

} // namespace duct


