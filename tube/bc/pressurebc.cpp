#include "pressurebc.h"
#include "conduction.h"


// Can be removed later
#include "energyeq.h"

namespace tube{
  dmat PressureBcEnergyEq::dUi(const TubeVertex& v) const {
    TRACE(40,"PressureBcEnergyEq::dUi()");
    dmat dUi=Energy::dUi(v);
    const Energy& e=static_cast<const Energy&>(*v.eq[2]);
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;

    d gamma=e.gamma(v);
    d gamfac=gamma/(gamma-1.0);

    const variable::var& pL=static_cast<const LeftPressure&>(v).pL;
    // dUi+=fDFT*pL.diagt()*iDFT;

    // dUi+=-v.wL0*gamfac*fDFT*pL.diagt()*iDFT;
    
    return dUi;
  }
  dmat PressureBcEnergyEq::dUip1(const TubeVertex& v) const {
    TRACE(40,"PressureBcEnergyEq::dUip1()");
    dmat dUip1=Energy::dUip1(v);
    const Energy& e=static_cast<const Energy&>(*v.eq[2]);
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;

    d gamma=e.gamma(v);
    d gamfac=gamma/(gamma-1.0);

    const variable::var& pL=static_cast<const LeftPressure&>(v).pL;
    // dUip1+=-fDFT*pL.diagt()*iDFT;

    // dUip1+=-v.wL1*gamfac*fDFT*pL.diagt()*iDFT;
    
    return dUip1;
  }  
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
    TRACE(100,"TL(0):"<<TL);
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
    if(eq[2]->getType()!=EqType::Ise){
      TRACE(100,"Changing energy equation...");
      // eq[2]=&peq;
    }
    LeftPressure::updateW(thisseg);

  }

  void LeftPressure::updateW(const SegBase& thisseg)
  {
    // Change continuity equation for an open boundary
    TRACE(80,"LeftPressure::updateW()");
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
    eWgim1= 0;
    eWgi=   wRl-wL0;
    eWgip1= wRr-wL1;

    eWkinim1=0;
    eWkini=wRl/pow(lg.SfR,2)-wL0/pow(lg.SfL,2);    
    eWkini=wRr/pow(lg.SfR,2)-wL1/pow(lg.SfL,2);    

    // Deprecated
    eWjim1 = 0;
    eWji   =-wRl;
    eWjip1 =-wRr;

    eWji+=wL0;
    eWjip1+=wL1;

    d x0=lg.vxi;
    d x1=lg.vxip1;
    d denom=x0*(1-x0/x1);
    TRACE(100,"denom:"<<denom);
    d x0_ov_x1sq=pow(x0/x1,2);
    eWc1=0;
    eWc2=lg.SfL/denom;
    eWc3=lg.SfR/lg.dxp;
    // eWc3=-lg.SfL*pow(xp1/xp2,2)/denom +  lg.SfR/dxp;
    eWc4=-lg.SfR/lg.dxp -lg.SfL*x0_ov_x1sq/denom;


    // TODO Fill this further!

  }
  vd LeftPressure::msource() const{
    TRACE(5,"LeftPressure::msource()");
    vd msource(gc->Ns,fillwith::zeros);
    msource=-1.0*lg.SfL*pL();
    // This one should not yet be scaled. The scaling is done in the
    // error term after adding this source.
    // TRACE(-1,"msource:"<<msource);
    return msource;
  }
  vd LeftPressure::esource() const {
    TRACE(5,"LeftPressure::esource()");
    vd esource=TubeVertex::esource();
    const dmat& fDFT=gc->fDFT;
    // TRACE(100,"Stupid hack to test energy source");
    vd TLt=TL.tdata();
    d TL0=TL(0);
    TRACE(100,"TL(0):"<<TL0);
    const Energy& e=static_cast<const Energy&>(*eq[2]);
    d gamma=e.gamma(*this);
    d gamfac=gamma/(gamma-1.0);
    
    // TRACE(100,"TL:"<<TL());

    d xp2=lg.vxip1;
    d xp1=lg.vxi;

    vd T0=gc->T0*vd(gc->Ns,fillwith::ones);
    vd kappaL=e.kappaL(*this);
    // esource+=lg.SfL*num*fDFT*(kappaL%TLt)/denom;
    d x0=lg.vxi;
    d x1=lg.vxip1;
    d denom=x0*(1.0-x0/x1);
    // TRACE(100,"Denom:"<<denom);
    d x0_ov_x1sq=pow(x0/x1,2);
    esource+=-1.0*(1-x0_ov_x1sq)*lg.SfL*fDFT*(kappaL%TLt)/denom;
    // esource+=fDFT*(U.tdata()%pL.tdata());
    // TRACE(100,"wL0:"<<wL0);
    // TRACE(100,"WL1:"<<wL1);
    // esource+=-gamfac*fDFT*(pL.tdata()%(wL0*U.tdata()+wL1*right->U.tdata()));
    
    return esource;    
  }

} // namespace tube
