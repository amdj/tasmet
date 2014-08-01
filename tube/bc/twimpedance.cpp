// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "twimpedance.h"
#include "tubevertex.h"


namespace tube{
  TwImpedance::TwImpedance(us segnr):TubeBcVertex(segnr),eright(*this){
    TRACE(8,"TwImpedance constructor");
    // Change continuity equation for open boundary
  }
  TwImpedance::TwImpedance(const TwImpedance& other):TwImpedance(other.segNumber())
  {
    TRACE(8,"TwImpedance copy cc.");
  }
  TwImpedance& TwImpedance::operator=(const TwImpedance& o){
    TRACE(8,"TwImpedance copy assignment operator");
    setSegNumber(o.segNumber());
    return *this;
  }
  vd TwImpedanceEnergyEq::error(const TubeVertex& v) const {
    return Energy::error(v);
  }
  void TwImpedance::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"TwImpedance::initTubeVertex(), vertex "<< i <<".");
    TubeVertex::initTubeVertex(i,thisseg);
    eq.at(1)=&mright;
    // eq[2]=&eright;
    // eq[2]=&is;
    // mright.Init(*thisseg.gc);
    // is.Init(*thisseg.gc);
    // eright.Init(*thisseg.gc);
    updateW(thisseg);
  }
  
  void TwImpedance::updateW(const SegBase& thisseg){
    const Geom& geom=thisseg.geom;
    TRACE(8,"TwImpedance::updateW()");
    cWddt=lg.vVf;
    cWim1=wRNm2-wLl;
    cWi  =wRNm1-wLr;
    cWip1=0;

    // Change momentum eq for open boundary
    mWddt=lg.vVf/lg.vSf;	// VERY IMPORTANT: update this!! (Is still zero)
    mWuim1=-wLl/lg.SfL+wRNm2/lg.SfR;
    mWui	=-wLr/lg.SfL+wRNm1/lg.SfR;
    mWuip1=0;

    mWpim1=-lg.SfL*wLl;
    mWpi	=-lg.SfL*wLr+(lg.SfL-lg.SfR);
    mWpip1=0;
    
    eWgim1=-wLl+wRNm2;
    eWgi  =-wLr+wRNm1;
    eWgip1=0;
    
    eWjim1=wLl-wRNm2;
    eWji  =wLr-wRNm1;
    eWjip1=0;
    
    d xi=geom.vx(i);
    d xim1=geom.vx(i-1);	 
    d dxm=xi-xim1;


    eWc1=-lg.SfL/lg.dxm;
    eWc2= lg.SfL/lg.dxm;
    eWc3= lg.SfR/lg.xr;
    eWc4=0;

    // Last but not least: point momentum and energy eq to new equation!
    // eq[1]=&mright;
    // eq[2]=&eright;        
    
  }

  vd TwImpedance::esource() const {
    // Source term related to temperature boundary condition
    TRACE(100,"TwImpedance::esource()");
    const dmat& fDFT=gc->fDFT;
    const dmat& iDFT=gc->iDFT;
    vd esource(gc->Ns,fillwith::zeros);
    d p0=gc->p0;
    d T0=gc->T0;
    d c0=gc->gas.cm(T0);
    d gamma=gc->gas.gamma(T0);
    // vd p0t=p0*vd(gc->Ns,fillwith::ones);
    vd urt=(wRNm1*U.tdata()+wRNm2*left->U.tdata())/lg.SfR;
    vd prt=p0*pow(1.0+((gamma-1.0)/2.0)*urt/c0,2.0*gamma/(gamma-1.0));
    // vd prt=(p0+wRNm1*p.tdata()+wRNm2*left->p.tdata());
    // TRACE(100,"prt:"<<prt);
    vd Trt=T0*pow(prt/p0,(gamma-1.0)/gamma);		// Adiabatic compression/expansion
    // TRACE(100,"Trt:"<<Trt);
    vd T0t=gc->T0*vd(gc->Ns,fillwith::ones);
    vd kappaR=gc->gas.kappa(T0t);
    // vd kappaR=eright.kappaR();
    // TRACE(100,"kappaR:"<<kappaR);    
    esource+=-1.0*(lg.SfR/xhalf)*fDFT*(kappaR%Trt);
    return esource;  
  }

  TwImpedanceEnergyEq::TwImpedanceEnergyEq(TwImpedance& twimp):impedancevertex(twimp){
    TRACE(4,"TwImpedanceEnergyEq::TwImpedanceEnergyEq()");
  }
  dmat TwImpedanceEnergyEq::dUi(const TubeVertex& v) const{
    TRACE(4,"TwImpedanceEnergyEq::dUi()");
    dmat dUi=Energy::dUi(v);
    assert(v.left!=NULL);
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;
    d xhalf=impedancevertex.xhalf;
    d T0=v.gc->T0;
    d p0=v.gc->p0;
    d c0=v.gc->gas.cm(T0);
    d gamma=v.gc->gas.gamma(T0);
    d z0=p0*gamma/c0;		// linearized characteristic impedance

    vd urt=(v.wRNm1*v.U.tdata()+v.wRNm2*v.left->U.tdata())/v.lg.SfR;
    // Find adiabatic temperature rize
    vd prt=p0*pow(1.0+((gamma-1.0)/2.0)*urt/c0,2.0*gamma/(gamma-1.0));
    vd Trt=T0*pow(prt/p0,(gamma-1.0)/gamma);

    vd kappaR=v.gc->gas.kappa(v.wRNm1*v.T.tdata()+v.wRNm2*v.left->T.tdata());
    vd Ztd=(z0/v.lg.SfR)*(pow(1.0+((gamma-1.0)/2.0)*urt/c0,(gamma+1.0)/(gamma-1.0))); // Impedance dp/dUr in time domain

    vd dTdprt=(T0/p0)*((gamma-1.0)/gamma)*pow(prt/p0,-1.0/gamma);
    dUi+=(-1.0*(v.lg.SfR/impedancevertex.xhalf)*v.wRNm1)*fDFT*diagmat(kappaR%dTdprt%Ztd)*iDFT;

    // TRACE(100,"dUi returns:"<<dUi);
    return dUi;
  }
  dmat TwImpedanceEnergyEq::dUim1(const TubeVertex& v) const {
    TRACE(4,"TwImpedanceEnergyEq::dUim1()");
    dmat dUim1=Energy::dUim1(v);
    assert(v.left!=NULL);
    dUim1+=TwImpedanceEnergyEq::dUi(v)*v.wRNm2/v.wRNm1;
    return dUim1;
  }

  TwImpedanceMomentumEq::TwImpedanceMomentumEq(){
    TRACE(6,"TwImpedanceMomentumEq::TwImpedanceMomentumEq()");
      }
  vd TwImpedanceMomentumEq::error(const TubeVertex& v) const {
    TRACE(20,"TwImpedanceMomentumEq::Error()");
    vd error(v.gc->Ns,fillwith::zeros);    
    // Add the normal stuff
    error+=Momentum::error(v);

    d T0=v.gc->T0;
    d c0=v.gc->gas.cm(T0);
    d p0=v.gc->p0;
    d gamma=v.gc->gas.gamma(T0);
    vd ur=(v.wRNm1*v.U()+v.wRNm2*v.left->U())/v.lg.SfR;
    vd urtd=v.gc->iDFT*ur;
    vd pr=v.gc->fDFT*(p0*pow(1.0+((gamma-1.0)/2.0)*urtd/c0,2.0*gamma/(gamma-1.0))-p0);
    vd errorZ=v.lg.SfR*pr;

    error+=errorZ;
    return error;
  }
  dmat TwImpedanceMomentumEq::dUi(const TubeVertex& v) const{
    TRACE(50,"TwImpedanceMomentumEq::dUi()");
    dmat dUi=Momentum::dUi(v);
    const us& Ns=v.gc->Ns;
    d T0=v.gc->T0;
    d c0=v.gc->gas.cm(T0);
    d p0=v.gc->p0;
    d gamma=v.gc->gas.gamma(T0);
    d z0=p0*gamma/c0;
    // TRACE(30,"z0:"<<z0);
    vd ur=(v.wRNm1*v.U()+v.wRNm2*v.left->U())/v.lg.SfR;
    vd urtd=v.gc->iDFT*ur;
    dmat Z=v.gc->fDFT*(z0/v.lg.SfR)*diagmat(pow(1.0+((gamma-1.0)/2.0)*urtd/c0,(gamma+1.0)/(gamma-1.0)))*v.gc->iDFT;

    dUi+=v.wRNm1*v.lg.SfR*Z;
    return dUi;
  }

  dmat TwImpedanceMomentumEq::dUim1(const TubeVertex& v) const{
    TRACE(50,"TwImpedanceMomentumEq::dUim1()");
    dmat dUim1=Momentum::dUim1(v);    
    const us& Ns=v.gc->Ns;
    d T0=v.gc->T0;
    d p0=v.gc->p0;
    d c0=v.gc->gas.cm(T0);
    d gamma=v.gc->gas.gamma(T0);
    d z0=p0*gamma/c0;

    vd ur=(v.wRNm1*v.U()+v.wRNm2*v.left->U())/v.lg.SfR;
    vd urtd=v.gc->iDFT*ur;
    dmat Z=v.gc->fDFT*(z0/v.lg.SfR)*diagmat(pow(1.0+((gamma-1.0)/2.0)*urtd/c0,(gamma+1.0)/(gamma-1.0)))*v.gc->iDFT;
    TRACE(-2,"Z:"<<Z);
    dUim1+=v.wRNm2*v.lg.SfR*Z;
    return dUim1;
  }

} // namespace tube












