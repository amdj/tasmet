// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
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
  void TwImpedance::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"TwImpedance::initTubeVertex(), vertex "<< i <<".");
    pr=var(gc);
    // pr.set(0);
    TubeVertex::initTubeVertex(i,thisseg);
    eqs.push_back(std::unique_ptr<TubeEquation>(twright.copy()));
    eqs.at(3).reset(midstate.copy()); // Replace perfect gas equation
    vars.push_back(&pr);
    TwImpedance::updateW(thisseg);
  }
  
  void TwImpedance::updateW(const SegBase& thisseg){
    TRACE(8,"TwImpedance::updateW()");
    cWim1=w.wRNm2-w.wLl;
    cWi  =w.wRNm1-w.wLr;
    cWip1=0;

    // Change momentum eq for open boundary
    mWddt=lg.vVf/lg.vSf;	// VERY IMPORTANT: update this!! (Is still zero)
    mWuim1=-w.wLl/lg.SfL+w.wRNm2/lg.SfR;
    mWui	=-w.wLr/lg.SfL+w.wRNm1/lg.SfR;
    mWuip1=0;
    WARN("Not yet updated!");
    // mWpim1=-lg.SfL*w.wLl;
    // mWpi	=-lg.SfL*w.wLr+(lg.SfL-lg.SfR);
    // mWpip1=0;
    
    eWgim1=-w.wLl;
    eWgim =-w.wLr;

    eWgUim1pR=w.wRNm2;
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


  vd RightTwImpedanceEq::error(const TubeVertex& v) const {
    TRACE(10,"TwImpedanceMomentumEq::Error()");
    vd error(v.gc->Ns,fillwith::zeros);    
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
    
    JacRow jac(dofnr,3);
    TRACE(15,"Dofnr:"<<dofnr);
    jac+=dUi(v);
    jac+=dUim1(v);
    jac+=dpR(v);

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
    JacCol dpR(v.pR());
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

} // namespace tube












