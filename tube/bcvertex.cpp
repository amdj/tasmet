// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "bcvertex.h"



namespace tube{
  LeftPressure::LeftPressure(const Tube& t):TubeVertex(t,0),pL(vop),TL(vop){
    pL(0)=tube.gc.p0;
    TL(0)=tube.gc.T0;
    Init();
  }
  LeftPressure::LeftPressure(const Tube& t,variable::var& pres,variable::var& temp):TubeVertex(t,0),pL(pres),TL(temp) {
    TRACE(0,"LeftPressure full constructor");
    Init();
  }
  LeftPressure::LeftPressure(const Tube& t,variable::var& pres,d T0):TubeVertex(t,0),pL(pres),TL(vop){
    TRACE(0,"LeftPressure constructor for given pressure. Temperature computed");    
    
    d gamma=e.gamma();
    d p0=pL(0);

    vd TLt=T0*pow((pL.tdata()/p0),(gamma-1.0)/gamma);		// Adiabatic compression/expansion
    TRACE(-1,"TLt:"<<pL.tdata());
    TL.settdata(TLt);
    Init();
  }
  void LeftPressure::Init(){
    // Change continuity equation for an open boundary
    c.Wi=c.wRl-c.wL0;
    c.Wip1=c.wRr-c.wL1;
    
    // Change momentum equation for open boundary, and prescribed pressure
    m.Wuim1=0;
    m.Wui=m.wRl-m.wL0;
    m.Wuip1=m.wRr-m.wL1;
    m.Wpi=m.wRl*m.SfR+m.SfL-m.SfR;
    m.Wpip1=m.SfR*m.wRr;
    // Change energy equation for open boundary and prescribed pressure
    e.Whim1=0;
    e.Whi=e.wRl-e.wL0;
    e.Whip1=e.wRr-e.wL1;

    e.Wjim1=0;
    e.Wji=e.wL0-e.wRl;
    e.Wjip1=e.wL1-e.wRr;

    d vxi=tube.geom.vx(0);
    d vxip1=tube.geom.vx(1);
    d dxp=vxip1-vxi;
    
    e.Wc1=0;
    e.Wc2=e.SfL/vxi;
    e.Wc3=e.SfR/dxp;
    e.Wc4=-e.SfR/dxp;
      // TODO Fill this further!

  }
  vd LeftPressure::msource() const{
    TRACE(0,"LeftPressure::msource()");
    vd msource(Ns,fillwith::ones);
    msource=-1.0*e.SfL*pL();
    return msource;
  }
  vd LeftPressure::esource() const {
    TRACE(0,"LeftPressure::esource()");
    vd esource(Ns,fillwith::zeros);
    d vxi=tube.geom.vx(0);
    vd TLt=TL.tdata();
    vd kappaL=tube.gas.kappa(TLt);
    // TRACE(10,"Important: put esource on when going back to full energy eq!");
    esource+=-1.0*e.SfL*(e.fDFT*(kappaL%TLt))/vxi;
    TRACE(-1,"esource:"<<esource);
    return esource;    
  }
  LeftPressure::~LeftPressure(){
    TRACE(-5,"LeftPressure destructor");
  }



} // namespace tube
