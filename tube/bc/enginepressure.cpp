#define TRACERPLUS (10)
#include "enginepressure.h"


namespace tube{


  // LeftEnginePressureContinuityEq& LeftEnginePressureContinuityEq::operator=(const Continuity& c){
  //   TRACE(15,"LeftEnginePressureContinuityEq::operator=()");
  //   Wddt=c.Wddt;
  //   Wim1=c.Wim1;
  //   Wi=c.Wi;
  //   Wip1=c.Wip1;
  //   Wart1=c.Wart1;
  //   Wart2=c.Wart2;
  //   Wart3=c.Wart3;
  //   Wart4=c.Wart4;
  // }
  
  vd LeftEnginePressureState::error(const TubeVertex& v) const {
    TRACE(15,"LeftEnginePressureState::error()");
    vd error(v.gc->Ns,fillwith::zeros);
    if(v.gc->Nf>0){
      error(1)=v.pL()(1)-lep->amplitude();
      error(2)=v.pL()(2);
    }
    return error;
  }
  JacRow LeftEnginePressureState::jac(const TubeVertex& v) const{
    TRACE(15,"LeftEnginePressureState::jac()");
    JacRow jac(dofnr,5);
    jac+=dpL(v);
    jac+=dTi(v);
    jac+=drhoi(v);
    if(v.left){
      jac+=dTim1(v);
      jac+=drhoim1(v);
    }
    else{
      jac+=dTip1(v);
      jac+=drhoip1(v);
    }

    return jac;
  }  

  JacCol LeftEnginePressureState::dpL(const TubeVertex& v)
   const {

    return JacCol(v.pL(),arma::eye(v.gc->Ns,v.gc->Ns));
  }

  JacCol LeftEnginePressureState::dTi(const TubeVertex& v)
   const {
    TRACE(0,"LeftEnginePressureState::dTi()");
    return  JacCol(v.T,v.zero);
  }
  JacCol LeftEnginePressureState::drhoi(const TubeVertex& v)
   const {
    TRACE(0,"LeftEnginePressureState::drhoi()");
    return JacCol(v.rho,v.zero);
  }
  JacCol LeftEnginePressureState::dTim1(const TubeVertex& v)
   const {
    TRACE(0,"LeftEnginePressureState::dTim1()");
    return JacCol(v.left->T,v.zero);
  }
  JacCol LeftEnginePressureState::drhoim1(const TubeVertex& v)
   const {
    TRACE(0,"LeftEnginePressureState::drhoim1()");
    return JacCol(v.left->rho,v.zero);
  }
  JacCol LeftEnginePressureState::dTip1(const TubeVertex& v)
   const {
    TRACE(0,"LeftEnginePressureState::dTip1()");
    return JacCol(v.right->T,v.zero);
  }
  JacCol LeftEnginePressureState::drhoip1(const TubeVertex& v)
   const {
    TRACE(0,"LeftEnginePressureState::drhoip1()");
    return JacCol(v.right->rho,v.zero);
  }
  
  
  
  void LeftEnginePressure::initTubeVertex(us i,const Tube& thistube){
    TRACE(15,"LeftEnginePressure::initTubeVertex()");

    lepc.init(thistube);
    lepm.init(thistube);
    lepe.init(thistube);
    leps.init(thistube);

    cA.init(thistube);
    mA.init(thistube);
    eA.init(thistube);
    sA.init(thistube);

    
    LeftAdiabaticWall::initTubeVertex(i,thistube);
    eqs.at(0)=&lepc;
    eqs.at(1)=&lepm;
    eqs.at(2)=&lepe;
    eqs.at(3)=&leps;

    cA.Wddt=c.Wddt;
    cA.Wim1=0;
    cA.Wi=wRl-wL0;
    cA.Wip1=wRr-wL1;

    const LocalGeom& rlg=right->lg;
    // Change momentum equation for open boundary, and prescribed pressure
    mA.Wddt=m.Wddt;
    mA.Wuim1=0;
    mA.Wui=wRl/lg.vSf-wL0/lg.SfL;
    mA.Wuip1=wRr/rlg.vSf-wL1/lg.SfL;
    mA.WpL=m.WpL;
    mA.WpR=m.WpR;

    eA.Wddt=e.Wddt;
    eA.Wddtkin=e.Wddtkin;    
    // Change energy equation for open boundary and prescribed pressure
    eA.Wgim1= 0;
    // Left boundary, velocity 
    eA.Wgim=-wL0;
    eA.WgUip1pL=-wL1;
    
    eA.Wgip=wRl;
    eA.Wgip1= wRr;    

    eA.Wkinim1=0;
    eA.Wkini=wRl/pow(lg.vSf,2)-wL0/pow(lg.SfL,2);    
    eA.Wkinip1=wRr/pow(rlg.vSf,2)-wL1/pow(lg.SfL,2);    
    
  }


  
}                // namespace tube
