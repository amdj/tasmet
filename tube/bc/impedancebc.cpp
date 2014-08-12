// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "impedancebc.h"
#include "tubevertex.h"
#include "w.h"

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
    updateW(thisseg);
  }
  
  void RightImpedance::updateW(const Tube& thisseg){
    TRACE(8,"RightImpedance::updateW()");
    cWddt=lg.vVf;
    cWim1=w.wRNm2-w.wLl;
    cWi  =w.wRNm1-w.wLr;
    cWip1=0;

    // Change momentum eq for open boundary
    mWddt=lg.vVf/lg.vSf;	// VERY IMPORTANT: update this!! (Is still zero)
    mWuim1=-w.wLl/w.vSfL+w.wRNm2/lg.SfR;
    mWui	=-w.wLr/w.vSfL+w.wRNm1/lg.SfR;
    mWuip1=0;

    mWpim1=-w.vSfL*w.wLl;
    mWpi	=-w.vSfL*w.wLr+(w.vSfL-lg.SfR);
    mWpip1=0;
    
    eWgim1=-w.wLl+w.wRNm2;
    eWgi  =-w.wLr+w.wRNm1;
    eWgip1=0;
    WARN("Not updated for kinetic energy terms!");
    eWc1=-w.vSfL/lg.dxm;
    eWc2= w.vSfL/lg.dxm;
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
    vd errorZ=v.lg.SfR*Z%(v.w.wRNm1*v.U()+v.w.wRNm2*v.left->U());
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
    dUi+=v.w.wRNm1*v.lg.SfR*diagmat(Z);
    // For pressure boundary condition
    // dUi.row(0).zeros();
    return dUi;
  }

  dmat RightImpedanceMomentumEq::dUim1(const TubeVertex& v) const {
    TRACE(1,"RightImpedanceMomentumEq::dUim1()");
    dmat dUim1=Momentum::dUim1(v);    
    dUim1+=v.w.wRNm2*v.lg.SfR*diagmat(Z);
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












