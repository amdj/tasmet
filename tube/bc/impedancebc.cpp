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

  void RightImpedance::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"RightImpedance::Init(), vertex "<< i <<".");
    TubeVertex::initTubeVertex(i,thisseg);
    eq.at(1)=&mright;
    updateW();
  }
  
  void RightImpedance::updateW(){
    TRACE(8,"RightImpedance::updateW()");
    cWddt=lg.vVf;
    cWim1=wRNm2-wLl;
    cWi  =wRNm1-wLr;
    cWip1=0;

    // Change momentum eq for open boundary
    mWddt=lg.vVf/lg.vSf;	// VERY IMPORTANT: update this!! (Is still zero)
    mWuim1=-wLl/lg.SfL+wRNm2/lg.SfR;
    mWui	=-wLr/lg.SfL+wRNm1/lg.SfR;
    mWuip1=0;

    mWpim1=-lg.SfL*wLl;
    mWpi	=-lg.SfL*wLr+(lg.SfL-lg.SfR);
    mWpip1=0;
    
    eWgim1=-wLl+wRNm2;
    eWgi  =-wLr+wRNm1;
    eWgip1=0;
    
    eWc1=-lg.SfL/lg.dxm;
    eWc2= lg.SfL/lg.dxm;
    eWc3=0;
    eWc4=0;

    // Conduction terms are not changed.
  }
  RightImpedanceMomentumEq::RightImpedanceMomentumEq(TubeBcVertex& tv,vd& Z):Z(Z){
    TRACE(6,"RightImpedanceMomentumEq::RightImpedanceMomentumEq()");
    TRACE(6,"Z:\n"<<Z)
      }
  vd RightImpedanceMomentumEq::error(const TubeVertex& v) const{
    TRACE(40,"RightImpedanceMomentumEq::Error()");
    vd error=Momentum::error(v);
    // SfR*p = SfR*Z*U
    vd errorZ=v.lg.SfR*Z%(v.wRNm1*v.U()+v.wRNm2*v.left->U());
    error+=errorZ;
    return error;
  }
  // dmat RightImpedanceMomentumEq::dpi(){
  //   TRACE(40,"RightImpedanceMomentumEq::dpi()");
  //   dmat dpi=Momentum::dpi();
  //   return dpi;
  // }
  dmat RightImpedanceMomentumEq::dUi(const TubeVertex& v) const {
    TRACE(40,"RightImpedanceMomentumEq::dUi()");
    dmat dUi=Momentum::dUi(v);
    dUi+=v.wRNm1*v.lg.SfR*diagmat(Z);
    // For pressure boundary condition
    // dUi.row(0).zeros();
    
    // For velocity boundary condition
    // dUi.row(0).zeros();
    // dUi(0,0)=MOM_SCALE0*MOM_SCALE;
      
    return dUi;
  }

  dmat RightImpedanceMomentumEq::dUim1(const TubeVertex& v) const {
    TRACE(1,"RightImpedanceMomentumEq::dUim1()");
    dmat dUim1=Momentum::dUim1(v);    
    dUim1+=v.wRNm2*v.lg.SfR*diagmat(Z);
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












