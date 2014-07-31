// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "impedancebc.h"
#include "tubevertex.h"
#include "momscale.h"


namespace tube{
  RightImpedance::RightImpedance(us segnr,vd Z1):TubeBcVertex(segnr),Z(Z1),mright(*this,Z){
    TRACE(8,"RightImpedance constructor");
    // Change continuity equation for open boundary
  }
  RightImpedance::RightImpedance(const RightImpedance& other):RightImpedance(other.segNumber(),other.Z)
  {
    TRACE(8,"RightImpedance copy cc."); 
  }
  RightImpedance& RightImpedance::operator=(const RightImpedance& o){
    TRACE(8,"RightImpedance copy assignment operator");
    setSegNumber(o.segNumber());
    Z=o.Z;
    return *this;
  }

  void RightImpedance::Init(us i,const SegBase& thisseg)
  {
    TRACE(8,"RightImpedance::Init(), vertex "<< i <<".");
    TubeVertex::Init(i,thisseg);
    eq[1]=&mright;
    mright.Init(*thisseg.gc);
    updateW(thisseg.geom);
  }
  
  void RightImpedance::updateW(const Geom& geom){
    TRACE(8,"RightImpedance::updateW()");
    c.Wddt=vVf;
    c.Wim1=wRNm2-wLl;
    c.Wi  =wRNm1-wLr;
    c.Wip1=0;

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

    // Conduction terms are not changed.
  }
  RightImpedanceMomentumEq::RightImpedanceMomentumEq(TubeBcVertex& tv,vd& Z):Momentum(tv),Z(Z){
    TRACE(6,"RightImpedanceMomentumEq::RightImpedanceMomentumEq()");
    TRACE(6,"Z:\n"<<Z)
      }
  vd RightImpedanceMomentumEq::Error(){
    TRACE(40,"RightImpedanceMomentumEq::Error()");
    vd error=Momentum::Error();
    // SfR*p = SfR*Z*U
    vd errorZ=v.SfR*Z%(v.wRNm1*v.U()+v.wRNm2*v.left->U());
    error+=errorZ;
    return error;
  }
  // dmat RightImpedanceMomentumEq::dpi(){
  //   TRACE(40,"RightImpedanceMomentumEq::dpi()");
  //   dmat dpi=Momentum::dpi();
  //   return dpi;
  // }
  dmat RightImpedanceMomentumEq::dUi(){
    TRACE(40,"RightImpedanceMomentumEq::dUi()");
    dmat dUi=Momentum::dUi();
    dUi+=v.wRNm1*v.SfR*diagmat(Z);
    // For pressure boundary condition
    // dUi.row(0).zeros();
    
    // For velocity boundary condition
    // dUi.row(0).zeros();
    // dUi(0,0)=MOM_SCALE0*MOM_SCALE;
      
    return dUi;
  }

  dmat RightImpedanceMomentumEq::dUim1(){
    TRACE(1,"RightImpedanceMomentumEq::dUim1()");
    dmat dUim1=Momentum::dUim1();    
    dUim1+=v.wRNm2*v.SfR*diagmat(Z);
    return dUim1;
  }

  // dmat RightImpedanceMomentumEq::dpim1(){
  //   TRACE(1,"RightImpedanceMomentumEq::dpim1()");
  //   dmat dpim1=Momentum::dpim1();
  //   // For velocity and pressure boundary condition
  //   // dpim1.row(0).zeros();
  //   return dpim1;
  // }
  // dmat RightImpedanceMomentumEq::drhoi(){
  //   TRACE(1,"RightImpedanceMomentumEq::drhoi()");
  //   dmat drhoi=Momentum::drhoi();
  //   // For velocity and pressure boundary condition
  //   // drhoi.row(0).zeros();
  //   return drhoi;
  // }
  // dmat RightImpedanceMomentumEq::drhoim1(){
  //   TRACE(1,"RightImpedanceMomentumEq::drhoim1()");
  //   dmat drhoim1=Momentum::drhoim1();
  //   // For velocity and pressure boundary condition
  //   // drhoim1.row(0).zeros();
  //   return drhoim1;
  // }

} // namespace tube












