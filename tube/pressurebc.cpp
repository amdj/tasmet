#include "pressurebc.h"
#include "bcvertex.h"
namespace tube{

  LeftPressure::LeftPressure(const Tube& t,variable::var& pres,variable::var& temp):TubeBcVertex(t,0),pL(pres),TL(temp) {
     TRACE(0,"LeftPressure full constructor");
   }
   LeftPressure::LeftPressure(const Tube& t,variable::var& pres):TubeBcVertex(t,0),pL(pres),TL(gc){
     TRACE(0,"LeftPressure constructor for given pressure. Temperature computed");    
     d T0=tube.gc.T0;
     d gamma=tube.gas.gamma(T0);
     vd p0(Ns,fillwith::ones); p0*=tube.gc.p0;
     // TRACE(-1,"p0:"<<p0);
     vd TLt=T0*pow((p0+pL.tdata())/p0,gamma/(gamma-1.0));		// Adiabatic compression/expansion
     // TRACE(-1,"TLt:"<<TLt);
     TL.settdata(TLt);

   }
   void LeftPressure::updateW(){
     // Change continuity equation for an open boundary
     TRACE(1,"LeftPressure::updateW()");
     TubeVertex::updateW();
     c.Wim1=0;
     c.Wi=wRl-wL0;
     c.Wip1=wRr-wL1;

     // Change momentum equation for open boundary, and prescribed pressure
     m.Wuim1=0;
     m.Wui=wRl-wL0;
     m.Wuip1=wRr-wL1;

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

     d vxi=tube.geom.vx(0);
     d vxip1=tube.geom.vx(1);
     d dxp=vxip1-vxi;

     e.Wc1=0;
     e.Wc2=SfL/vxi;
     e.Wc3=SfR/dxp;
     e.Wc4=-SfR/dxp;
       // TODO Fill this further!

   }
   vd LeftPressure::msource() const{
     TRACE(0,"LeftPressure::msource()");
     vd msource(Ns,fillwith::ones);
     msource=-1.0*SfL*pL();
     // This one should not yet be scaled. The scaling is done in the
     // error term after adding this source.
    TRACE(-1,"msource:"<<msource);
    return msource;
  }
  vd LeftPressure::esource() const {
    TRACE(0,"LeftPressure::esource()");
    vd esource(Ns,fillwith::zeros);
    d vxi=tube.geom.vx(0);
    vd TLt=TL.tdata();
    vd kappaL=tube.gas.kappa(TLt);
    // TRACE(10,"Important: put esource on when going back to full energy eq!");
    esource+=-1.0*SfL*(e.fDFT*(kappaL%TLt))/vxi;
    TRACE(-1,"esource:"<<esource);
    return esource;    
  }

} // namespace tube
