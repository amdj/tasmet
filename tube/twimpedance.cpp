// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "twimpedance.h"
#include "tubevertex.h"
#include "momscale.h"


namespace tube{
  TwImpedance::TwImpedance(us segnr):TubeBcVertex(segnr),mright(*this){
    TRACE(8,"TwImpedance constructor");
    eq[1]=&mright;
    // Change continuity equation for open boundary
  }
  TwImpedance::TwImpedance(const TwImpedance& other):TwImpedance(other.segNumber())
  {
    TRACE(8,"TwImpedance copy cc."); 
  }
  TwImpedance& TwImpedance::operator=(const TwImpedance& o){
    TRACE(8,"TwImpedance copy assignment operator");
    setSegNumber(o.segNumber());
    eq[1]=&mright;
    return *this;
  }

  void TwImpedance::Init(us i,const SegBase& thisseg)
  {
    TRACE(8,"TwImpedance::Init(), vertex "<< i <<".");
    TubeVertex::Init(i,thisseg);
    eq[1]=&mright;
    mright.Init(*thisseg.gc);
    updateW(thisseg);
  }
  
  void TwImpedance::updateW(const SegBase& thisseg){
    const Geom& geom=thisseg.geom;
    TRACE(8,"TwImpedance::updateW()");
    c.Wddt=vVf;
    c.Wim1=wRNm2-wLl;
    c.Wi  =wRNm1-wLr;
    c.Wip1=0;
    // TRACE(20,"SfL:"<<SfL);
    // TRACE(20,"SfR:"<<SfR);
    // TRACE(20,"wRNm1:"<<wRNm1);
    // TRACE(20,"wRNm2:"<<wRNm2);

    // Change momentum eq for open boundary
    mright.Wddt=vVf/vSf;	// VERY IMPORTANT: update this!! (Is still zero)
    mright.Wuim1=-wLl/SfL+wRNm2/SfR;
    mright.Wui	=-wLr/SfL+wRNm1/SfR;
    mright.Wuip1=0;

    mright.Wpim1=-SfL*wLl;
    mright.Wpi	=-SfL*wLr+(SfL-SfR);
    mright.Wpip1=0;
    
    e.Wgim1=-wLl+wRNm2;
    e.Wgi  =-wLr+wRNm1;
    e.Wgip1=0;
    
    e.Wjim1=wLl-wRNm2;
    e.Wji  =wLr-wRNm1;
    e.Wjip1=0;
    
    d xi=geom.vx(i);
    d xim1=geom.vx(i-1);	 
    d dxm=xi-xim1;

    xhalf=xR-vxi;	       
    e.Wc1=-SfL/dxm;
    e.Wc2= SfL/dxm;
    e.Wc3=SfR/xhalf;
    e.Wc4=0;

    // Last but not least: point momentum eq to new equation!
    eq[1]=&mright;
    // Conduction terms are not changed.
    
  }
  vd TwImpedance::esource(){
    // Source term related to temperature boundary condition
    TRACE(6,"TwImpedance::esource()");
    const dmat& fDFT=gc->fDFT;
    vd esource(gc->Ns,fillwith::zeros);
    variable::var pR=wRNm2*left->p+wRNm1*p;
    
    d T0=gc->T0;
    d gamma=gc->gas.gamma(T0);
    vd p0(gc->Ns,fillwith::ones); p0*=gc->p0;
    // TRACE(-1,"p0:"<<p0);
    vd TRt=T0*pow((p0+pR.tdata())/p0,gamma/(gamma-1.0));		// Adiabatic compression/expansion
    vd kappaR=gc->gas.kappa(TRt);

    esource+=-1.0*SfR*fDFT*(kappaR%TRt)/xhalf;
    return esource;  
  }
  
  TwImpedanceMomentumEq::TwImpedanceMomentumEq(TubeBcVertex& tv):Momentum(tv){
    TRACE(6,"TwImpedanceMomentumEq::TwImpedanceMomentumEq()");
      }
  vd TwImpedanceMomentumEq::Error(){
    TRACE(4,"TwImpedanceMomentumEq::Error()");

    vd error(v.gc->Ns,fillwith::zeros);    
    // Add the normal stuff
    error+=Momentum::Error();

    d T0=v.gc->T0;
    d c0=v.gc->gas.cm(T0);
    d p0=v.gc->p0;
    d gamma=v.gc->gas.gamma(T0);
    vd ur=(v.wRNm1*v.U()+v.wRNm2*v.left->U())/v.SfR;
    vd urtd=v.gc->iDFT*ur;
    vd pr=v.gc->fDFT*(p0*pow(1.0+((gamma-1.0)/2.0)*urtd/c0,2.0*gamma/(gamma-1.0))-p0);
    // TRACE(100,"pr:"<<pr);
    vd errorZ=MOM_SCALE*v.SfR*pr;

    errorZ(0)*=MOM_SCALE0;
    error+=errorZ;
    return error;
  }
  dmat TwImpedanceMomentumEq::dpi(){
    TRACE(1,"TwImpedanceMomentumEq::dpi()");
    dmat dpi=Momentum::dpi();
    return dpi;
  }
  dmat TwImpedanceMomentumEq::dUi(){
    TRACE(1,"TwImpedanceMomentumEq::dUi()");
    dmat dUi=Momentum::dUi();
    const us& Ns=v.gc->Ns;
    d T0=v.gc->T0;
    d c0=v.gc->gas.cm(T0);
    d p0=v.gc->p0;
    d gamma=v.gc->gas.gamma(T0);
    d z0=p0*gamma/c0;
    // TRACE(30,"z0:"<<z0);
    vd ur=(v.wRNm1*v.U()+v.wRNm2*v.left->U())/v.SfR;
    vd urtd=v.gc->iDFT*ur;
    
    dmat Z=v.gc->fDFT*(z0/v.SfR)*diagmat(pow(1.0+((gamma-1.0)/2.0)*urtd/c0,(gamma+1.0)/(gamma-1.0)))*v.gc->iDFT;
    dmat adddUi=MOM_SCALE*v.wRNm1*v.SfR*Z;
    adddUi.row(0)*=MOM_SCALE0;
    dUi+=adddUi;
    return dUi;
  }

  dmat TwImpedanceMomentumEq::dUim1(){
    TRACE(1,"TwImpedanceMomentumEq::dUim1()");
    dmat dUim1=Momentum::dUim1();    
    const us& Ns=v.gc->Ns;
    d T0=v.gc->T0;
    d p0=v.gc->p0;
    d c0=v.gc->gas.cm(T0);
    d gamma=v.gc->gas.gamma(T0);
    d z0=p0*gamma/c0;
    // TRACE(30,"z0:"<<z0);

    vd ur=(v.wRNm1*v.U()+v.wRNm2*v.left->U())/v.SfR;
    vd urtd=v.gc->iDFT*ur;
    
    dmat Z=v.gc->fDFT*(z0/v.SfR)*diagmat(pow(1.0+((gamma-1.0)/2.0)*urtd/c0,(gamma+1.0)/(gamma-1.0)))*v.gc->iDFT;
    dmat adddUim1=MOM_SCALE*v.wRNm2*v.SfR*Z;
    adddUim1.row(0)*=MOM_SCALE0;
    dUim1+=adddUim1;
    return dUim1;
  }

  dmat TwImpedanceMomentumEq::dpim1(){
    TRACE(1,"TwImpedanceMomentumEq::dpim1()");
    dmat dpim1=Momentum::dpim1();
    return dpim1;
  }
  dmat TwImpedanceMomentumEq::drhoi(){
    TRACE(1,"TwImpedanceMomentumEq::drhoi()");
    dmat drhoi=Momentum::drhoi();
    // For velocity and pressure boundary condition
    // drhoi.row(0).zeros();
    return drhoi;
  }
  dmat TwImpedanceMomentumEq::drhoim1(){
    TRACE(1,"TwImpedanceMomentumEq::drhoim1()");
    dmat drhoim1=Momentum::drhoim1();
    // For velocity and pressure boundary condition
    // drhoim1.row(0).zeros();
    return drhoim1;
  }

} // namespace tube












