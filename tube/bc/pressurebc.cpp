#include "pressurebc.h"
#include "conduction.h"
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
    vd TLt=T0*pow((p0+pL.tdata())/p0,(gamma-1.0)/gamma);		// Adiabatic compression/expansion
    
    TL.settdata(TLt);
    TRACE(100,"TL(0):"<<TL(0));
  }
  LeftPressure::LeftPressure(const LeftPressure& other):LeftPressure(other.segNumber(),other.pL,other.TL)
  {
    TRACE(8,"LeftPressure copy constructor");
  }
  void LeftPressure::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"LeftPressure::Init()");
    pL.gc=thisseg.gc;
    TL.gc=thisseg.gc;
    TubeVertex::initTubeVertex(i,thisseg);
    LeftPressure::updateW(thisseg);
  }

  void LeftPressure::updateW(const SegBase& thisseg)
  {
    // Change continuity equation for an open boundary
    TRACE(8,"LeftPressure::updateW()");
    const Geom& geom=thisseg.geom;
    cWim1=0;
    cWi=wRl-wL0;
    cWip1=wRr-wL1;

    // Change momentum equation for open boundary, and prescribed pressure
    mWuim1=0;
    mWui=wRl/lg.SfR-wL0/lg.SfL;
    mWuip1=wRr/lg.SfR-wL1/lg.SfL;

    mWpim1=0;     
    mWpi=wRl*lg.SfR+lg.SfL-lg.SfR;
    mWpip1=lg.SfR*wRr;
    // Change energy equation for open boundary and prescribed pressure
    eWgim1=0;
    eWgi=wRl-wL0;
    eWgip1=wRr-wL1;

    eWjim1=0;
    eWji=wL0-wRl;
    eWjip1=wL1-wRr;

    d vxi=geom.vx(0);
    d vxip1=geom.vx(1);
    d dxp=vxip1-vxi;

    d xp2=vxip1;
    d xp1=vxi;

    d denom=xp1*(1-xp1/xp2);
    
    eWc1=0;
    eWc2=lg.SfL/denom;
    
    eWc3=-lg.SfL*pow(xp1/xp2,2)/denom +  lg.SfR/dxp;
    eWc4=-lg.SfR/dxp;

    
    // TODO Fill this further!

  }
  vd LeftPressure::msource() const{
    TRACE(2,"LeftPressure::msource()");
    vd msource(gc->Ns,fillwith::zeros);
    msource=-1.0*lg.SfL*pL();
    // This one should not yet be scaled. The scaling is done in the
    // error term after adding this source.
    TRACE(-1,"msource:"<<msource);
    return msource;
  }
  vd LeftPressure::esource() const {
    TRACE(100,"LeftPressure::esource()");
    vd esource(gc->Ns,fillwith::zeros);
    const dmat& fDFT=gc->fDFT;
    vd TLt=TL.tdata()+50;
    TRACE(100,"TL:"<<TL());

    d xp2=lg.vxip1;
    d xp1=lg.vxi;

    d denom=xp1*(1-xp1/xp2);
    
    d num=-(1-pow(xp1/xp2,2));    
    
    vd T0=gc->T0*vd(gc->Ns,fillwith::ones);
    vd kappaL=gc->gas.kappa(T0);
    esource+=lg.SfL*num*fDFT*(kappaL%TLt)/denom;
    return esource;    
  }

} // namespace tube
