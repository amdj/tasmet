// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "impedancebc.h"
#include "tubevertex.h"
#include "momscale.h"


namespace tube{
  RightImpedance::RightImpedance(us segnr,vd Z1):TubeBcVertex(segnr),Z(Z1),mright(*this,Z){
    TRACE(8,"RightImpedance constructor");
    eq[1]=&mright;
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
    eq[1]=&mright;
    return *this;
  }

  void RightImpedance::Init(us i,const Globalconf& gc,const Geom& geom)
  {
    TRACE(8,"RightImpedance::Init(), vertex "<< i <<".");
    TubeVertex::Init(i,gc,geom);
    eq[1]=&mright;
    mright.Init(gc);
    updateW(geom);
  }
  
  void RightImpedance::updateW(const Geom& geom){
    TRACE(8,"RightImpedance::updateW()");
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
  RightImpedanceMomentumEq::RightImpedanceMomentumEq(TubeBcVertex& tv,vd& Z):Momentum(tv),Z(Z){
    TRACE(6,"RightImpedanceMomentumEq::RightImpedanceMomentumEq()");
    TRACE(6,"Z:\n"<<Z)
      }
  vd RightImpedanceMomentumEq::Error(){
    TRACE(4,"RightImpedanceMomentumEq::Error()");
    vd error(gc->Ns,fillwith::zeros);
    // Add the normal stuff
    error+=Momentum::Error();  
    vd errorZ=MOM_SCALE*vertex.SfR*Z%(vertex.wRNm1*vertex.U()+vertex.wRNm2*vertex.left->U());

    errorZ(0)*=MOM_SCALE0;
    error+=errorZ;

    return error;
  }
  dmat RightImpedanceMomentumEq::dpi(){
    TRACE(1,"RightImpedanceMomentumEq::dpi()");
    dmat dpi=Momentum::dpi();
    return dpi;
  }
  dmat RightImpedanceMomentumEq::dUi(){
    TRACE(1,"RightImpedanceMomentumEq::dUi()");
    dmat dUi=Momentum::dUi();
    dmat adddUi=MOM_SCALE*vertex.wRNm1*vertex.SfR*diagmat(Z);
    adddUi.row(0)*=MOM_SCALE0;
    dUi+=adddUi;

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
    dmat adddUim1=MOM_SCALE*vertex.wRNm2*vertex.SfR*diagmat(Z);
    adddUim1.row(0)*=MOM_SCALE0;
    dUim1+=adddUim1;
    // For velocity boundary condition
    // dUim1.row(0).zeros();
    return dUim1;
  }

  dmat RightImpedanceMomentumEq::dpim1(){
    TRACE(1,"RightImpedanceMomentumEq::dpim1()");
    dmat dpim1=Momentum::dpim1();
    // For velocity and pressure boundary condition
    // dpim1.row(0).zeros();
    return dpim1;
  }
  dmat RightImpedanceMomentumEq::drhoi(){
    TRACE(1,"RightImpedanceMomentumEq::drhoi()");
    dmat drhoi=Momentum::drhoi();
    // For velocity and pressure boundary condition
    // drhoi.row(0).zeros();
    return drhoi;
  }
  dmat RightImpedanceMomentumEq::drhoim1(){
    TRACE(1,"RightImpedanceMomentumEq::drhoim1()");
    dmat drhoim1=Momentum::drhoim1();
    // For velocity and pressure boundary condition
    // drhoim1.row(0).zeros();
    return drhoim1;
  }

} // namespace tube












