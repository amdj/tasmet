// file: bccell.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "impedancebc.h"
#include "cell.h"


namespace tube{
  RightImpedance::RightImpedance(vd Z1):TubeBcCell(),Z(Z1),mright(*this,Z){
    TRACE(8,"RightImpedance constructor");
    // Change continuity equation for open boundary
  }
  RightImpedance::RightImpedance(const RightImpedance& other):RightImpedance(other.Z)
  {
    TRACE(8,"RightImpedance copy cc."); 
  }
  RightImpedance& RightImpedance::operator=(const RightImpedance& o){
    TRACE(8,"RightImpedance copy assignment operator");
    Z=o.Z;
    return *this;
  }

  void RightImpedance::initCell(us i,const Tube& thisseg)
  {
    TRACE(8,"RightImpedance::Init(), cell "<< i <<".");
    Cell::initCell(i,thisseg);
    eqs.at(1)=&mright;
    eqs.at(1)->init(thisseg);
    updateW(thisseg);
  }
  
  void RightImpedance::updateW(const Tube& thisseg){
    TRACE(8,"RightImpedance::updateW()");
    c.Wddt=lg.vVf;
    c.Wim1=wRNm2-wLl;
    c.Wi  =wRNm1-wLr;
    c.Wip1=0;

    // Change momentum eq for open boundary
    m.Wddt=lg.vVf/lg.vSf;	// VERY IMPORTANT: update this!! (Is still zero)
    m.Wuim1=-wLl/vSfL+wRNm2/lg.SfR;
    m.Wui	=-wLr/vSfL+wRNm1/lg.SfR;
    m.Wuip1=0;
    WARN("Not yet updated!");
    // mWpim1=-vSfL*wLl;
    // mWpi	=-vSfL*wLr+(vSfL-lg.SfR);
    // mWpip1=0;
    
    e.Wgim1=-wLl;
    e.Wgim =-wLr;

    e.WgUim1pR=wRNm2;
    e.Wgip=wRNm1;

    e.Wgip1=0;
    
    WARN("Not updated for kinetic energy terms!");
    e.Wc1=-vSfL/lg.dxm;
    e.Wc2= vSfL/lg.dxm;
    e.Wc3=0;
    e.Wc4=0;

    // Conduction terms are not changed.
  }
  RightImpedanceMomentumEq::RightImpedanceMomentumEq(TubeBcCell& tv,vd& Z):Z(Z){
    TRACE(6,"RightImpedanceMomentumEq::RightImpedanceMomentumEq()");
    TRACE(6,"Z:\n"<<Z)
      }
  vd RightImpedanceMomentumEq::error(const Cell& v) const{
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
  JacCol RightImpedanceMomentumEq::dUi(const Cell& v) const {
    TRACE(40,"RightImpedanceMomentumEq::dUi()");
    JacCol dUi=Momentum::dUi(v);
    dUi+=v.wRNm1*v.lg.SfR*diagmat(Z);
    // For pressure boundary condition
    // dUi.row(0).zeros();
    return dUi;
  }

  JacCol RightImpedanceMomentumEq::dUim1(const Cell& v) const {
    TRACE(1,"RightImpedanceMomentumEq::dUim1()");
    JacCol dUim1=Momentum::dUim1(v);    
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












