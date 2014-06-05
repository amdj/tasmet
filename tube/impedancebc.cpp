// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "impedancebc.h"
#include "tubevertex.h"
#include "momscale.h"


namespace tube{
  RightImpedance::RightImpedance(const Tube& t,vd Z1):TubeBcVertex(t,t.getNcells()-1),Z(Z1),mright(t,*this,Z){
    TRACE(5,"RightImpedance constructor");
    // Change continuity equation for open boundary
  }
  void RightImpedance::updateW(){
    TRACE(1,"RightImpedance::updateW()");
    TubeVertex::updateW();
    c.Wddt=vVf;
    c.Wim1=wRNm2-wLl;
    c.Wi  =wRNm1-wLr;
    c.Wip1=0;
    // TRACE(20,"SfL:"<<SfL);
    // TRACE(20,"SfR:"<<SfR);
    // TRACE(20,"wRNm1:"<<wRNm1);
    // TRACE(20,"wRNm2:"<<wRNm2);

    // Change momentum eq for open boundary
    mright.Wddt=vVf/vSf;
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
    
    d xi=tube.geom.vx(i);
    d xim1=tube.geom.vx(i-1);	 
    d dxm=xi-xim1;
    e.Wc1=-SfL/dxm;
    e.Wc2= SfL/dxm;
    e.Wc3=0;
    e.Wc4=0;

    // Last but not least: point momentum eq to new equation!
    eq[1]=&mright;
    // Conduction terms are not changed.
    
  }
  RightImpedanceMomentumEq::RightImpedanceMomentumEq(const Tube& t,TubeBcVertex& tv,vd& Z):Momentum(t,tv),Z(Z){
    TRACE(20,"RightImpedanceMomentumEq::RightImpedanceMomentumEq()");
    TRACE(20,"Z:\n"<<Z)
  }
  vd RightImpedanceMomentumEq::Error(){
    TRACE(20,"RightImpedanceMomentumEq::Error()");
    vd error(Ns,fillwith::zeros);
    // Add the normal stuff
    error+=Momentum::Error();  
    // And add the right pressure as being an impedance times the
    // velocity:
    // TRACE(10,"wRNm1:"<<vertex.wRNm1);
    // TRACE(10,"wRNm2:"<<vertex.wRNm2);
    // TRACE(10,"Left U:"<<vertex.left->U());
    // TRACE(10,"This U:"<<vertex.U());    
    vd errorZ=MOM_SCALE*vertex.SfR*Z%(vertex.wRNm1*vertex.U()+vertex.wRNm2*vertex.left->U());
    // TRACE(10,"ErrorZ:\n"<<errorZ);
    errorZ(0)*=MOM_SCALE0;
    error+=errorZ;
    // TRACE(10,"Z:"<<Z);
    // TRACE(10,"Error:\n"<<error);
    // Use time-averaged momentum equation to fix the pressure to p0
    // at the last node
    // Impedance-like
    // error(0)=1.0*MOM_SCALE0*MOM_SCALE*(p(0)-Z(0)*U(0));

    // Just set velocity at last node equalt to zero :Velocity zero
    // error(0)=1.0*MOM_SCALE0*MOM_SCALE*U(0);

    // Pressure opening
    // error(0)+=wRNm1*MOM_SCALE0SfR*p(0)+wRNm2*SfR*left->p(0);

    
    // Pressure zero
    // error(0)=MOM_SCALE0*MOM_SCALE*p(0);


    return error;
  }
  dmat RightImpedanceMomentumEq::dpi(){
    TRACE(10,"RightImpedanceMomentumEq::dpi()");
    dmat dpi=Momentum::dpi();
    // dpi.row(0).zeros();

    // Set mean pressure to zero on last node
    // dpi(0,0)=MOM_SCALE*MOM_SCALE0;
    return dpi;
  }
  dmat RightImpedanceMomentumEq::dUi(){
    TRACE(10,"RightImpedanceMomentumEq::dUi()");
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
    TRACE(10,"RightImpedanceMomentumEq::dUim1()");
    dmat dUim1=Momentum::dUim1();    
    dmat adddUim1=MOM_SCALE*vertex.wRNm2*vertex.SfR*diagmat(Z);
    adddUim1.row(0)*=MOM_SCALE0;
    dUim1+=adddUim1;
    // For velocity boundary condition
    // dUim1.row(0).zeros();
    return dUim1;
  }

  dmat RightImpedanceMomentumEq::dpim1(){
    TRACE(10,"RightImpedanceMomentumEq::dpim1()");
    dmat dpim1=Momentum::dpim1();
    // For velocity and pressure boundary condition
    // dpim1.row(0).zeros();
    return dpim1;
  }
  dmat RightImpedanceMomentumEq::drhoi(){
    TRACE(10,"RightImpedanceMomentumEq::drhoi()");
    dmat drhoi=Momentum::drhoi();
    // For velocity and pressure boundary condition
    // drhoi.row(0).zeros();
    return drhoi;
  }
  dmat RightImpedanceMomentumEq::drhoim1(){
    TRACE(10,"RightImpedanceMomentumEq::drhoim1()");
    dmat drhoim1=Momentum::drhoim1();
    // For velocity and pressure boundary condition
    // drhoim1.row(0).zeros();
    return drhoim1;
  }

} // namespace tube












