// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
// #define TRACERPLUS 20

#include "tube.h"
#include "twimpedance.h"
#include "tubevertex.h"
#include "w.h"

namespace tube{
  TwImpedance::TwImpedance():TubeBcVertex(){
    TRACE(8,"TwImpedance constructor");
    // Change continuity equation for open boundary
  }
  TwImpedance::TwImpedance(const TwImpedance& other):TwImpedance()
  {
    TRACE(8,"TwImpedance copy cc.");
  }
  TwImpedance& TwImpedance::operator=(const TwImpedance& o){
    TRACE(8,"TwImpedance copy assignment operator");
    return *this;
  }
  // vd TwImpedanceEnergyEq::error(const TubeVertex& v) const {
  //   return Energy::error(v);
  // }
  void TwImpedance::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(12,"TwImpedance::initTubeVertex(), vertex "<< i <<".");
    TubeVertex::initTubeVertex(i,thisseg);
    // mright.init(thisseg);
    // eqs.at(1).reset(mright.copy());
    pr=var(gc);
    // sr.init(thisseg);
    // s.init(thisseg);
    righttwimp.init(thisseg);
    // eqs.at(3).reset(s.copy());
    eqs.push_back(std::unique_ptr<TubeEquation>(righttwimp.copy()));
    // eqs.push_back(std::unique_ptr<TubeEquation>(sr.copy()));
    vars.push_back(&pr);

    TwImpedance::updateW(thisseg);
  }
  
  void TwImpedance::updateW(const SegBase& thisseg){
    TRACE(12,"TwImpedance::updateW()");
    cWim1=w.wRNm2-w.wLl;
    cWi  =w.wRNm1-w.wLr;
    cWip1=0;

    // Change momentum eq for open boundary
    mWddt=lg.vVf/lg.vSf;	// VERY IMPORTANT: update this!! (Is still zero)
    mWuim1=-w.wLl/lg.SfL+w.wRNm2/lg.SfR;
    mWui	=-w.wLr/lg.SfL+w.wRNm1/lg.SfR;
    mWuip1=0;

    
    eWgim1=-w.wLl;
    eWgUim1pR=w.wRNm2;
    eWgim  =-w.wLr;
    eWgip=w.wRNm1;
    eWgip1=0;

    d SfLsq=pow(w.vSfL,2);
    d SfRsq=pow(lg.SfR,2);
    eWkinim1=-0.5*w.wLl/SfLsq+0.5*w.wRNm2/SfRsq;
    eWkini=-0.5*w.wLr/SfLsq+0.5*w.wRNm1/SfRsq;
    eWkinip1=0;

    eWc1=-w.vSfL/w.dxm;
    eWc2= lg.vSf/w.dxm;
    eWc3=lg.SfR/lg.xr;
    eWc4=0;
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
    // vd urt=(wRNm1*U.tdata()+wRNm2*left->U.tdata())/lg.SfR;
    // vd prt=p0*pow(1.0+((gamma-1.0)/2.0)*urt/c0,2.0*gamma/(gamma-1.0));
    vd prt=(p0+w.wRNm1*p.tdata()+w.wRNm2*left->p.tdata());
    // TRACE(100,"prt:"<<prt);
    vd Trt=T0*pow(prt/p0,(gamma-1.0)/gamma);		// Adiabatic compression/expansion
    // TRACE(100,"Trt:"<<Trt);
    vd T0t=gc->T0*vd(gc->Ns,fillwith::ones);
    vd kappaR=gc->gas.kappa(T0t);
    // vd kappaR=eright.kappaR();
    // TRACE(100,"kappaR:"<<kappaR);

    esource+=-1.0*(lg.SfR/lg.xr)*fDFT*(kappaR%Trt);
    return esource;  
  }

  // vd TwImpedanceMomentumEq::error(const TubeVertex& v) const {
  //   TRACE(10,"TwImpedanceMomentumEq::Error()");
  //   vd error(v.gc->Ns,fillwith::zeros);    
  //   // Add the normal stuff
  //   error+=Momentum::error(v);
  //   error-=v.mWpR*v.pR()();
  //   d T0=v.gc->T0;
  //   d c0=v.gc->gas.cm(T0);
  //   d p0=v.gc->p0;
  //   d gamma=v.gc->gas.gamma(T0);
  //   vd ur=(v.w.wRNm1*v.U()+v.w.wRNm2*v.left->U())/v.lg.SfR;
  //   vd urtd=v.gc->iDFT*ur;
  //   vd pr=v.gc->fDFT*(p0*pow(1.0+((gamma-1.0)/2.0)*urtd/c0,2.0*gamma/(gamma-1.0))-p0);
  //   vd errorZ=v.lg.SfR*pr;

  //   error+=errorZ;
  //   return error;
  // }
  // JacCol TwImpedanceMomentumEq::dpR(const TubeVertex& v) const{
  //   TRACE(10,"TwImpedanceMomentumEq::dpR() ----- not adding anythin");
  //   JacCol dpR(v.pR());
  //   dpR.setToAdd(false);
  //   return dpR;
  // }
  // JacCol TwImpedanceMomentumEq::dUi(const TubeVertex& v) const{
  //   TRACE(10,"TwImpedanceMomentumEq::dUi()");
  //   JacCol dUi=Momentum::dUi(v);
  //   const us& Ns=v.gc->Ns;
  //   d T0=v.gc->T0;
  //   d c0=v.gc->gas.cm(T0);
  //   d p0=v.gc->p0;
  //   d gamma=v.gc->gas.gamma(T0);
  //   d z0=p0*gamma/c0;
  //   // TRACE(30,"z0:"<<z0);
  //   vd ur=(v.w.wRNm1*v.U()+v.w.wRNm2*v.left->U())/v.lg.SfR;
  //   vd urtd=v.gc->iDFT*ur;
  //   dmat Z=v.gc->fDFT*(z0/v.lg.SfR)*diagmat(pow(1.0+((gamma-1.0)/2.0)*urtd/c0,(gamma+1.0)/(gamma-1.0)))*v.gc->iDFT;

  //   dUi+=v.w.wRNm1*v.lg.SfR*Z;
  //   return dUi;
  // }

  // JacCol TwImpedanceMomentumEq::dUim1(const TubeVertex& v) const{
  //   TRACE(10,"TwImpedanceMomentumEq::dUim1()");
  //   JacCol dUim1=Momentum::dUim1(v);    
  //   const us& Ns=v.gc->Ns;
  //   d T0=v.gc->T0;
  //   d p0=v.gc->p0;
  //   d c0=v.gc->gas.cm(T0);
  //   d gamma=v.gc->gas.gamma(T0);
  //   d z0=p0*gamma/c0;

  //   vd ur=(v.w.wRNm1*v.U()+v.w.wRNm2*v.left->U())/v.lg.SfR;
  //   vd urtd=v.gc->iDFT*ur;
  //   dmat Z=v.gc->fDFT*(z0/v.lg.SfR)*diagmat(pow(1.0+((gamma-1.0)/2.0)*urtd/c0,(gamma+1.0)/(gamma-1.0)))*v.gc->iDFT;
  //   TRACE(-2,"Z:"<<Z);
  //   dUim1+=v.w.wRNm2*v.lg.SfR*Z;
  //   return dUim1;
  // }
  vd RightTwImpedanceEq::error(const TubeVertex& v) const {
    TRACE(10,"TwImpedanceMomentumEq::Error()");
    vd error=StateR::error(v);
    // Add the normal stuff
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;    

    d T0=v.gc->T0;
    d c0=v.gc->gas.cm(T0);
    d p0=v.gc->p0;
    d rho0=v.gc->rho0;
    d gamma=v.gc->gas.gamma(T0);
    d z0=rho0*c0;

    error+=(v.w.wRNm1*v.U()+v.w.wRNm2*v.left->U())/v.lg.SfR; // Adding u

    d powfacp=(gamma-1.0)/(2*gamma);
    d prefac=-2*c0/(gamma-1.0);
    error+=prefac*fDFT*(pow((v.pR().tdata()+p0)/p0,powfacp)-1.0);
    return error;
  }
  JacRow RightTwImpedanceEq::jac(const TubeVertex& v) const{
    TRACE(10,"RightTwImpedanceEq::jac()");
    
    JacRow jac(dofnr,8);
    TRACE(15,"Dofnr:"<<dofnr);
    jac+=dUi(v);
    jac+=dUim1(v);
    jac+=dpR(v);
    jac+=drhoi(v);
    jac+=dTi(v);    
    jac+=drhoim1(v);
    jac+=dTim1(v);
    return jac;
  }
  JacCol RightTwImpedanceEq::dUi(const TubeVertex& v) const{
    TRACE(10,"RightTwImpedanceEq::dUi()");
    TRACE(50,"dUi wRNm1:"<< v.w.wRNm1);

    return JacCol(v.U,v.w.wRNm1*eye(v.gc->Ns,v.gc->Ns));
  }
  JacCol RightTwImpedanceEq::dUim1(const TubeVertex& v) const{
    TRACE(10,"RightTwImpedanceEq::dUim1()");
    return JacCol(v.left->U,v.w.wRNm2*eye(v.gc->Ns,v.gc->Ns));
  }
  JacCol RightTwImpedanceEq::dpR(const TubeVertex& v) const{
    TRACE(10,"RightTwImpedanceEq::dpR()");
    JacCol dpR=StateR::dpR(v);
    const us& Ns=v.gc->Ns;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;    
    d T0=v.gc->T0;
    d c0=v.gc->gas.cm(T0);
    d p0=v.gc->p0;
    d gamma=v.gc->gas.gamma(T0);
    d z0=p0*gamma/c0;
    // TRACE(30,"z0:"<<z0);
    d powfac=(gamma-2.0)/(2*gamma);
    dpR+=(-1/z0)*fDFT*diagmat(pow((v.pR().tdata()+p0)/p0,powfac))*iDFT;
    return dpR;
  }


  // TwImpedanceEnergyEq::TwImpedanceEnergyEq(TwImpedance& twimp):impedancevertex(twimp){
  //   TRACE(4,"TwImpedanceEnergyEq::TwImpedanceEnergyEq()");
  // }
  // dmat TwImpedanceEnergyEq::dUi(const TubeVertex& v) const{
  //   TRACE(4,"TwImpedanceEnergyEq::dUi()");
  //   dmat dUi=Energy::dUi(v);
  //   assert(v.left!=NULL);
  //   const dmat& fDFT=v.gc->fDFT;
  //   const dmat& iDFT=v.gc->iDFT;
  //   d xhalf=impedancevertex.xhalf;
  //   d T0=v.gc->T0;
  //   d p0=v.gc->p0;
  //   d c0=v.gc->gas.cm(T0);
  //   d gamma=v.gc->gas.gamma(T0);
  //   d z0=p0*gamma/c0;		// linearized characteristic impedance

  //   vd urt=(v.wRNm1*v.U.tdata()+v.wRNm2*v.left->U.tdata())/v.lg.SfR;
  //   // Find adiabatic temperature rize
  //   vd prt=p0*pow(1.0+((gamma-1.0)/2.0)*urt/c0,2.0*gamma/(gamma-1.0));
  //   vd Trt=T0*pow(prt/p0,(gamma-1.0)/gamma);

  //   vd kappaR=v.gc->gas.kappa(v.wRNm1*v.T.tdata()+v.wRNm2*v.left->T.tdata());
  //   vd Ztd=(z0/v.lg.SfR)*(pow(1.0+((gamma-1.0)/2.0)*urt/c0,(gamma+1.0)/(gamma-1.0))); // Impedance dp/dUr in time domain

  //   vd dTdprt=(T0/p0)*((gamma-1.0)/gamma)*pow(prt/p0,-1.0/gamma);
  //   dUi+=(-1.0*(v.lg.SfR/impedancevertex.xhalf)*v.wRNm1)*fDFT*diagmat(kappaR%dTdprt%Ztd)*iDFT;

  //   // TRACE(100,"dUi returns:"<<dUi);
  //   return dUi;
  // }
  // dmat TwImpedanceEnergyEq::dUim1(const TubeVertex& v) const {
  //   TRACE(4,"TwImpedanceEnergyEq::dUim1()");
  //   dmat dUim1=Energy::dUim1(v);
  //   assert(v.left!=NULL);
  //   dUim1+=TwImpedanceEnergyEq::dUi(v)*v.wRNm2/v.wRNm1;
  //   return dUim1;
  // }

} // namespace tube












