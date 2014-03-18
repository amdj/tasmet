#include "continuityeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{

  vd Continuity::Error(){	// Current error in continuity equation
    vd error(Ns,fillwith::zeros);
    error+=vVf*DDTfd*vertex.rho();
    vd rhoip1=tube.vvertex[i+1].rho.tdata();
    vd Uip1=tube.vvertex[i+1].U.tdata();
    vd rhoim1=tube.vvertex[i-1].rho.tdata();
    vd Uim1=tube.vvertex[i-1].U.tdata();
    error+=wRr*fDFT*(rhoip1%Uip1);
    error+=-1.0*wLl*fDFT*(rhoim1%Uim1);
    error+=(wRl-wLr)*fDFT*(vertex.rho.tdata()%vertex.U.tdata());
    return error;
  }
  dmat Continuity::drhoi(){
    TRACE(0,"Continuity::drhoi()");
    dmat drhoi=vVf*DDTfd;		// Initialize and add first term
    drhoi+=(wRl-wLr)*fDFT*diagtmat(vertex.U)*iDFT;
    return drhoi;
  }
  dmat Continuity::dUi(){
    return (wRl-wLr)*fDFT*diagtmat(vertex.rho)*iDFT;
  }
  dmat Continuity::dUip1(){
    TRACE(0,"Continuity::dUip1()");
    return wRr*fDFT*diagtmat(tube.vvertex[i+1].rho)*iDFT;
  }
  dmat Continuity::dUim1(){
    TRACE(0,"Continuity::dUim1()");
    return wLl*fDFT*diagtmat(tube.vvertex[i-1].rho)*iDFT;
  }
  dmat Continuity::drhoip1(){
    TRACE(0,"Continuity::drhoip1()");
    return wRr*fDFT*diagtmat(tube.vvertex[i+1].U)*iDFT;
  }
  dmat Continuity::drhoim1(){
    TRACE(0,"Continuity::drhoim1()");
    return wLl*fDFT*diagtmat(tube.vvertex[i-1].U)*iDFT;
  }

  Continuity::~Continuity(){}

} // Namespace tube
