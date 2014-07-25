#include "pressurebc.h"

namespace tube{

  LeftPressure::LeftPressure(us segnr,const var& pres,const var& temp):TubeBcVertex(segnr),pL(pres),TL(temp) {
    TRACE(8,"LeftPressure full constructor");
  }
  LeftPressure::LeftPressure(us segnr,const var& pres):TubeBcVertex(segnr),pL(pres),TL(*pres.gc){
    TRACE(8,"LeftPressure constructor for given pressure. Temperature computed");    
    const Globalconf* gc=pres.gc;
    d T0=gc->T0;
    d gamma=gc->gas.gamma(T0);
    vd p0(gc->Ns,fillwith::ones); p0*=gc->p0;
    // TRACE(-1,"p0:"<<p0);
    vd TLt=T0*pow((p0+pL.tdata())/p0,gamma/(gamma-1.0));		// Adiabatic compression/expansion
    // TRACE(-1,"TLt:"<<TLt);
    TL.settdata(TLt);

  }
  LeftPressure::LeftPressure(const LeftPressure& other):LeftPressure(other.segNumber(),other.pL,other.TL)
  {
    TRACE(8,"LeftPressure copy constructor");
    TRACE(8,"pL:"<<pL);
    TRACE(9,"LeftPressure left pointer:"<<left);
  }
  void LeftPressure::Init(us i,const SegBase& thisseg)
  {
    TRACE(8,"LeftPressure::Init()");
    pL.gc=thisseg.gc;
    TL.gc=thisseg.gc;
    TubeVertex::Init(i,thisseg);
    LeftPressure::updateW(thisseg);
  }

  void LeftPressure::updateW(const SegBase& thisseg)
  {
    // Change continuity equation for an open boundary
    TRACE(8,"LeftPressure::updateW()");
    const Geom& geom=thisseg.geom;
    c.Wim1=0;
    c.Wi=wRl-wL0;
    c.Wip1=wRr-wL1;

    // Change momentum equation for open boundary, and prescribed pressure
    m.Wuim1=0;
    m.Wui=wRl/SfR-wL0/SfL;
    m.Wuip1=wRr/SfR-wL1/SfL;

    m.Wpim1=0;     
    m.Wpi=wRl*SfR+SfL-SfR;
    m.Wpip1=SfR*wRr;
    // Change energy equation for open boundary and prescribed pressure
    e.Wgim1=0;
    e.Wgi=wRl-wL0;
    e.Wgip1=wRr-wL1;

    e.Wjim1=0;
    e.Wji=wL0-wRl;
    e.Wjip1=wL1-wRr;

    d vxi=geom.vx(0);
    d vxip1=geom.vx(1);
    d dxp=vxip1-vxi;

    e.Wc1=0;
    e.Wc2=SfL/vxi;
    e.Wc3=SfR/dxp;
    e.Wc4=-SfR/dxp;
    // TODO Fill this further!

  }
  vd LeftPressure::msource() const{
    TRACE(2,"LeftPressure::msource()");
    vd msource(gc->Ns,fillwith::zeros);
    msource=-1.0*SfL*pL();
    // This one should not yet be scaled. The scaling is done in the
    // error term after adding this source.
    TRACE(-1,"msource:"<<msource);
    return msource;
  }
  vd LeftPressure::esource() const {
    TRACE(2,"LeftPressure::esource()");
    vd esource(gc->Ns,fillwith::zeros);
    vd TLt=TL.tdata();
    vd kappaL=gc->gas.kappa(TLt);
    // TRACE(10,"Important: put esource on when going back to full energy eq!");
    esource+=-1.0*SfL*(gc->fDFT*(kappaL%TLt))/vxi;
    TRACE(-1,"esource:"<<esource);
    return esource;    
  }

} // namespace tube
