#include "tube.h"
#include "pressurebc.h"


// Can be removed later
#include "energyeq.h"
namespace tube{
  
  // dmat PressureBcEnergyEq::dUi(const TubeVertex& v) const {
  //   TRACE(40,"PressureBcEnergyEq::dUi()");
  //   dmat dUi=Energy::dUi(v);
  //   const Energy& e=static_cast<const Energy&>(*v.eq[2]);
  //   const dmat& DDTfd=v.gc->DDTfd;
  //   const dmat& fDFT=v.gc->fDFT;
  //   const dmat& iDFT=v.gc->iDFT;

  //   d gamma=e.gamma(v);
  //   d gamfac=gamma/(gamma-1.0);

  //   const variable::var& pL=static_cast<const LeftPressure&>(v).pL;
  //   // dUi+=fDFT*pL.diagt()*iDFT;

  //   // dUi+=-v.wL0*gamfac*fDFT*pL.diagt()*iDFT;
    
  //   return dUi;
  // }
  // dmat PressureBcEnergyEq::dUip1(const TubeVertex& v) const {
  //   TRACE(40,"PressureBcEnergyEq::dUip1()");
  //   dmat dUip1=Energy::dUip1(v);
  //   const Energy& e=static_cast<const Energy&>(*v.eq[2]);
  //   const dmat& DDTfd=v.gc->DDTfd;
  //   const dmat& fDFT=v.gc->fDFT;
  //   const dmat& iDFT=v.gc->iDFT;

  //   d gamma=e.gamma(v);
  //   d gamfac=gamma/(gamma-1.0);

  //   const variable::var& pL=static_cast<const LeftPressure&>(v).pL;
  //   // dUip1+=-fDFT*pL.diagt()*iDFT;

  //   // dUip1+=-v.wL1*gamfac*fDFT*pL.diagt()*iDFT;
    
  //   return dUip1;
  // }  

  LeftPressure::LeftPressure(const var& pres,const var& temp):TubeBcVertex(),pL(pres),TL(temp) {
    TRACE(8,"LeftPressure full constructor");
  }
  LeftPressure::LeftPressure(const var& pres):TubeBcVertex(),pL(pres),TL(*pres.gc){
    TRACE(8,"LeftPressure constructor for given pressure. Temperature computed");    
    const Globalconf* gc=pres.gc;
    d T0=gc->T0;
    d gamma=gc->gas.gamma(T0);
    vd p0(gc->Ns,fillwith::ones); p0*=gc->p0;
    vd TLt=T0*pow((p0+pL.tdata())/p0,(gamma-1.0)/gamma);		// Adiabatic compression/expansion
    TL.settdata(TLt);
    // TRACE(100,"TL(0):"<<TL);
  }
  LeftPressure::LeftPressure(const LeftPressure& other):LeftPressure(other.pL,other.TL)
  {
    TRACE(8,"LeftPressure copy constructor");
  }
  void LeftPressure::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"LeftPressure::initTubeVertex()");
    TubeVertex::initTubeVertex(i,thisseg);
    pL.gc=thisseg.gc;
    TL.gc=thisseg.gc;
    // TRACE(100,"TL:"<<TL());
    if(eq.at(2)->getType()!=EqType::Ise){
      TRACE(100,"Changing energy equation...");
      // eq[2]=&peq;
    }
    LeftPressure::updateW(thisseg);

  }
  void LeftPressure::show() const {
    cout << "LeftPressure boundary condition set at left side of tube.\n";
    TubeVertex::show();
  }
  void LeftPressure::updateW(const SegBase& thisseg)
  {
    // Change continuity equation for an open boundary
    TRACE(8,"LeftPressure::updateW()");

    cWim1=0;
    cWi=w.wRl-w.wL0;
    cWip1=w.wRr-w.wL1;

    // Change momentum equation for open boundary, and prescribed pressure
    mWuim1=0;
    mWui=w.wRl/w.vSfR-w.wL0/w.vSfL;
    mWuip1=w.wRr/w.vSfR-w.wL1/w.vSfL;

    mWpim1=0;     
    mWpi=w.wRl*w.vSfR+w.vSfL-w.vSfR;
    mWpip1=w.vSfR*w.wRr;
    // Change energy equation for open boundary and prescribed pressure
    eWgim1= 0;
    eWgi=   w.wRl-w.wL0;
    eWgip1= w.wRr-w.wL1;

    eWkinim1=0;
    eWkini=w.wRl/pow(w.vSfR,2)-w.wL0/pow(w.vSfL,2);    
    eWkinip1=w.wRr/pow(w.vSfR,2)-w.wL1/pow(w.vSfL,2);    

    TRACE(100,"Sf:"<<w.vSf);
    d x0=w.xvi;
    d x1=w.xvip1;
    TRACE(100,"x0:"<<x1);
    TRACE(100,"x1:"<<x1);
    d denom=x0*(1-x0/x1);
    TRACE(100,"denom:"<<denom);
    d x0_ov_x1sq=pow(x0/x1,2);
    eWc1=0;
    // eWc2=lg.SfL/denom;
    eWc2=lg.SfL/x0;
    eWc3=lg.vSf/w.dxp;
    eWc4=-w.vSfR/w.dxp;    
    // eWc3=-w.vSfL*pow(x0/x1,2)/denom +  lg.vSf/w.dxp;
    // eWc4=-w.vSfR/w.dxp -lg.SfL*x0_ov_x1sq/denom;

  }
  vd LeftPressure::msource() const{
    TRACE(5,"LeftPressure::msource()");
    vd msource(gc->Ns,fillwith::zeros);
    msource=-1.0*lg.SfL*pL();
    return msource;
  }
  vd LeftPressure::esource() const {
    TRACE(5,"LeftPressure::esource()");
    vd esource=TubeVertex::esource();
    const dmat& fDFT=gc->fDFT;
    // TRACE(100,"Stupid hack to test energy source");
    vd TLt=TL.tdata();
    d TL0=TL(0);
    // TRACE(100,"TL(0):"<<TL0);
    const Energy& e=static_cast<const Energy&>(*eq[2]);
    d gamma=e.gamma(*this);
    d gamfac=gamma/(gamma-1.0);
    
    // TRACE(100,"TL:"<<TL());
    const LocalGeom rlg(lg.geom->localGeom(i+1));
    d xp2=rlg.xvi;
    d xp1=lg.xvi;

    vd T0=gc->T0*vd(gc->Ns,fillwith::ones);
    vd kappaL=e.kappaL(*this);
    // esource+=lg.SfL*num*fDFT*(kappaL%TLt)/denom;
    d x0=lg.xvi;
    d x1=lg.xvi+w.dxp;
    d denom=x0*(1.0-x0/x1);
    // TRACE(100,"Denom:"<<denom);
    d x0_ov_x1sq=pow(x0/x1,2);
    TRACE(12,"esource:"<<esource);
    // esource+=-1.0*(1-x0_ov_x1sq)*lg.SfL*fDFT*(kappaL%TLt)/denom;
    esource+=-1.0*lg.SfL*fDFT*(kappaL%TLt)/x0;
    // TRACE(100,"esourcefac:"<<-1.0*lg.SfL/x0);
    // TRACE(100,"esourcefac:"<<-1.0*(1-x0_ov_x1sq)*lg.SfL/denom);
    TRACE(12,"esource:"<<esource);
    // esource+=fDFT*(U.tdata()%pL.tdata());
    // TRACE(100,"wL0:"<<wL0);
    // TRACE(100,"WL1:"<<wL1);
    // esource+=-gamfac*fDFT*(pL.tdata()%(wL0*U.tdata()+wL1*right->U.tdata()));
    
    return esource;    
  }

} // namespace tube



