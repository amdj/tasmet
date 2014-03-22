#include "energyeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{

  Energy::Energy(const Tube& tube,const TubeVertex& gp):Equation(tube,gp){
    TRACE(0,"Energy constructor done");
  }
  vd Energy::Error(){		// Error in momentum equation
    TRACE(0,"Energy::Error()");
    vd error(Ns,fillwith::zeros);

    d gamma=this->gamma();
    TRACE(-1,"gamma: " << gamma);


    vd ptip1=tube.vvertex[i+1].p.tdata();
    vd ptim1=tube.vvertex[i-1].p.tdata();

    vd Utip1=tube.vvertex[i+1].U.tdata();
    vd Utim1=tube.vvertex[i-1].U.tdata();
    vd Uti=vertex.U.tdata();
    vd pti=vertex.p.tdata();

    error+=vVf*DDTfd*vertex.p()/(gamma-1.0);
    error+=(wRl-wLr)*fDFT*(pti%Uti)/(gamma-1.0);
    error+=wRr*fDFT*(ptip1%((gamma/(gamma-1))*Utip1-Uti));
    error+=-1.0*wLl*fDFT*(ptim1%((gamma/(gamma-1.0))*Utim1-Uti));

    // right boundary
    d xip1=tube.geom.vx(i+1);
    d xim1=tube.geom.vx(i-1);
    d xi=tube.geom.vx(i);
    d dxp=xip1-xi;
    d dxm=xi-xim1;
    vd Tip1=tube.vvertex[i+1].T.tdata();
    vd Tim1=tube.vvertex[i-1].T.tdata();
    vd Ti=vertex.T.tdata();
    error+=     SfL*fDFT*((kappaL()/dxm)%(Ti-Tim1));
    error+=-1.0*SfR*fDFT*((kappaR()/dxp)%(Tip1-Ti));
    return error;
  }
  dmat Energy::dpi(){
    TRACE(0,"Energy::dpi()");
    d T0=vertex.T(0);
    d gamma=tube.gas.gamma(T0);
    dmat dpi=zero;
    dpi+=(vVf/(gamma-1.0))*DDTfd;
    dpi+=(gamma*(wRl-wLr)/(gamma-1.0))*fDFT*diagtmat(vertex.U)*iDFT;
    return dpi;
  }

  dmat Energy::dUi(){
    TRACE(0,"Energy::dUi()");
    dmat dUi=zero;			    // Initialize with zeros
    d T0=vertex.T(0);
    d gamma=tube.gas.gamma(T0);
    dUi+=(gamma*(wRl-wLr)/(gamma-1.0))*fDFT*diagtmat(vertex.p)*iDFT;
    dUi+=fDFT*diagmat(wLl*tube.vvertex[i-1].p.tdata()-wRr*tube.vvertex[i+1].p.tdata())*iDFT;
    return dUi;
  }
  dmat Energy::dTi(){
    TRACE(0,"Energy::dTi()");
    dmat dTi=zero;

    d xip1=tube.geom.vx(i+1);
    d xim1=tube.geom.vx(i-1);
    d xi=tube.geom.vx(i);
    d dxp=xip1-xi;
    d dxm=xi-xim1;
    
    dTi+=fDFT*((SfL/dxm)*diagmat(kappaL())+(SfR/dxp)*diagmat(kappaR()))*iDFT;
    return dTi;
  }
  dmat Energy::dpip1(){
    TRACE(0,"Energy::dpip1()");
    vd Utip1=tube.vvertex[i+1].U.tdata();
    vd Uti=vertex.U.tdata();
    dmat dpip1=zero;
    d gamma=this->gamma();
    dpip1+=wRr*fDFT*diagmat((gamma/(gamma-1))*Utip1-Uti)*iDFT;
    return dpip1;
  }
  dmat Energy::dUip1(){
    TRACE(0,"Energy::dUip1()");
    d gamma=this->gamma();
    vd ptip1=tube.vvertex[i+1].p.tdata();
    dmat dUip1=zero;
    dUip1+=(gamma/(gamma-1.0)*fDFT*diagmat(ptip1))*iDFT;
    return dUip1;
  }
  dmat Energy::dTip1(){
    TRACE(0,"Energy::dTip1()");
    dmat dTip1=zero;

    d xip1=tube.geom.vx(i+1);
    d xim1=tube.geom.vx(i-1);
    d xi=tube.geom.vx(i);
    d dxp=xip1-xi;
    d dxm=xi-xim1;
    
    dTip1+=-1.0*fDFT*diagmat(SfR*kappaR()/(dxp))*iDFT;
    return dTip1;
  }
  dmat Energy::dpim1(){
    TRACE(0,"Energy::dpim1()");
    vd Utim1=tube.vvertex[i-1].U.tdata();
    vd Uti=vertex.U.tdata();
    dmat dpip1=zero;
    d gamma=this->gamma();
    dpip1+=-wLl*fDFT*diagmat((gamma/(gamma-1))*Utim1-Uti)*iDFT;
    return dpip1;
  }
  dmat Energy::dUim1(){
    TRACE(0,"Energy::dUim1()");
    d gamma=this->gamma();
    vd ptim1=tube.vvertex[i-1].p.tdata();
    dmat dUim1=zero;
    dUim1+=-1.0*(gamma/(gamma-1.0)*fDFT*diagmat(ptim1))*iDFT;
    return dUim1;
    }
  dmat Energy::dTim1(){
    TRACE(0,"Energy::dTim1()");
    dmat dTim1=zero;

    d xip1=tube.geom.vx(i+1);
    d xim1=tube.geom.vx(i-1);
    d xi=tube.geom.vx(i);
    d dxp=xip1-xi;
    d dxm=xi-xim1;
    
    dTim1+=fDFT*diagmat(SfL*kappaL()/(dxm))*iDFT;
    return dTim1;
  }    
  vd Energy::kappaL(){
    vd Tti=vertex.T.tdata();
    vd kappait=tube.gas.kappa(Tti);
    vd Ttim1=tube.vvertex[i-1].T.tdata();
    vd kappaitm1=tube.gas.kappa(Ttim1);
    vd kappaL=wLr*kappait+wLl*kappaitm1;	// Conductivity at the left boundary
    return kappaL;
  }
  vd Energy::kappaR(){
    vd Tti=vertex.T.tdata();
    vd kappait=tube.gas.kappa(Tti);
    vd Ttip1=tube.vvertex[i+1].T.tdata();
    vd kappaitp1=tube.gas.kappa(Ttip1);    
    vd kappaR=wRl*kappait+wRr*kappaitp1;	
    return kappaR;
  }
  d Energy::gamma(){
    d T0=vertex.T(0);
    return tube.gas.gamma(T0);
  }
  Energy::~Energy(){
    TRACE(-5,"Energy destructor");
  }
    
} // namespace tube




