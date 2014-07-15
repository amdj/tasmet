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

  void TwImpedance::Init(us i,const Globalconf& gc,const Geom& geom)
  {
    TRACE(8,"TwImpedance::Init(), vertex "<< i <<".");
    TubeVertex::Init(i,gc,geom);
    eq[1]=&mright;
    mright.Init(gc);
    updateW(geom);
  }
  
  void TwImpedance::updateW(const Geom& geom){
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
    e.Wc1=-SfL/dxm;
    e.Wc2= SfL/dxm;
    e.Wc3=0;
    e.Wc4=0;

    // Last but not least: point momentum eq to new equation!
    eq[1]=&mright;
    // Conduction terms are not changed.
    
  }
  TwImpedanceMomentumEq::TwImpedanceMomentumEq(TubeBcVertex& tv):Momentum(tv){
    TRACE(6,"TwImpedanceMomentumEq::TwImpedanceMomentumEq()");
      }
  vd TwImpedanceMomentumEq::Error(){
    TRACE(4,"TwImpedanceMomentumEq::Error()");

    vd error(gc->Ns,fillwith::zeros);    
    // Add the normal stuff
    error+=Momentum::Error();

    d T0=gc->T0;
    d c0=gc->gas.cm(T0);
    d p0=gc->p0;
    d gamma=gc->gas.gamma(T0);
    vd ur=(vertex.wRNm1*vertex.U()+vertex.wRNm2*vertex.left->U())/vertex.SfR;
    vd urtd=gc->iDFT*ur;
    vd pr=gc->fDFT*(p0*pow(1.0+((gamma-1.0)/2.0)*urtd/c0,2.0*gamma/(gamma-1.0))-p0);
    // TRACE(100,"pr:"<<pr);
    vd errorZ=MOM_SCALE*vertex.SfR*pr;

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
    const us& Ns=gc->Ns;
    d T0=gc->T0;
    d c0=gc->gas.cm(T0);
    d p0=gc->p0;
    d gamma=gc->gas.gamma(T0);
    d z0=p0*gamma/c0;
    // TRACE(30,"z0:"<<z0);
    vd ur=(vertex.wRNm1*vertex.U()+vertex.wRNm2*vertex.left->U())/vertex.SfR;
    vd urtd=gc->iDFT*ur;
    
    dmat Z=gc->fDFT*(z0/vertex.SfR)*diagmat(pow(1.0+((gamma-1.0)/2.0)*urtd/c0,(gamma+1.0)/(gamma-1.0)))*gc->iDFT;
    dmat adddUi=MOM_SCALE*vertex.wRNm1*vertex.SfR*Z;
    adddUi.row(0)*=MOM_SCALE0;
    dUi+=adddUi;
    return dUi;
  }

  dmat TwImpedanceMomentumEq::dUim1(){
    TRACE(1,"TwImpedanceMomentumEq::dUim1()");
    dmat dUim1=Momentum::dUim1();    
    const us& Ns=gc->Ns;
    d T0=gc->T0;
    d p0=gc->p0;
    d c0=gc->gas.cm(T0);
    d gamma=gc->gas.gamma(T0);
    d z0=p0*gamma/c0;
    // TRACE(30,"z0:"<<z0);

    vd ur=(vertex.wRNm1*vertex.U()+vertex.wRNm2*vertex.left->U())/vertex.SfR;
    vd urtd=gc->iDFT*ur;
    
    dmat Z=gc->fDFT*(z0/vertex.SfR)*diagmat(pow(1.0+((gamma-1.0)/2.0)*urtd/c0,(gamma+1.0)/(gamma-1.0)))*gc->iDFT;
    dmat adddUim1=MOM_SCALE*vertex.wRNm2*vertex.SfR*Z;
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












