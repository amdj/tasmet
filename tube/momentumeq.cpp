#include "momentumeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{

  Momentum::Momentum(const Tube& tube,const TubeVertex& gp):Equation(tube,gp){
    TRACE(0,"Momentum constructor done");
  }
  dmat Momentum::operator()(){
    TRACE(0,"Momentum::operator()");
    return Equation::operator()();
  }
  vd Momentum::Error(){		// Error in momentum equation
    TRACE(0,"Momentum::Error()");
    vd error(Ns,fillwith::zeros);
    error+=vVf*DDTfd*fDFT*(vertex.U.tdata()%vertex.rho.tdata())/vSf;

    // WATCH IT! ABOVE TERMS ARE ALL TIME DOMAIN DATA!!
    vd rhoi=vertex.rho.tdata();
    vd Ui=vertex.U.tdata();
    vd rhoip1=tube.vvertex[i+1].rho.tdata();
    vd rhoim1=tube.vvertex[i-1].rho.tdata();
    vd Uip1=tube.vvertex[i+1].U.tdata();
    vd Uim1=tube.vvertex[i-1].U.tdata();
    
    error+=wRr*fDFT*(rhoip1%Uip1%Uip1)/vSf;
    error+=(wRl-wLr)*fDFT*(rhoi%Ui%Ui);
    vd presterm(Ns,fillwith::zeros);
    presterm+=wRr*tube.vvertex[i+1].p()*tube.geom.Sf(i+1);    
presterm+=
    
    // error+=-1.0*tube.vvertex.at(i-1).p()*tube.geom.Sf(i-1)/(dxp+dxm);
    error+=tube.drag(i);
    
    error+=fDFT*(rhotip1%Uip1%Uip1)/(tube.geom.Sf(i+1)*(dxp+dxm));
    // error+=-1.0*fDFT*(rhotim1%Uim1%Uim1)/(tube.geom.Sf(i-1)*(dxp+dxm));
    return error;
  }
  dmat Momentum::drhoi(){
    TRACE(0,"Momentum::drhoi()");
    dmat Utdiag(Ns,Ns,fillwith::zeros);
    Utdiag.diag()=tube.vvertex.at(i).U.tdata();
    return DDTfd*fDFT*Utdiag*iDFT;
  }
  dmat Momentum::dUi(){
    TRACE(0,"Momentum::dUi()");
    TRACE(-1,"Ns:"<<Ns);
    dmat rhotdiag(Ns,Ns,fillwith::zeros);
    rhotdiag.diag()=tube.vvertex.at(i).rho.tdata();
    TRACE(-1,"rhotdiag:"<<rhotdiag);
    dmat ddragdU=tube.drag.dUi(i);
    return ddragdU+DDTfd*fDFT*rhotdiag*iDFT;
  }
  dmat Momentum::drhoip1(){
    // Todo: add this term!
    return zero;
  }
  dmat Momentum::dUip1(){
    // Todo: add this term!
    return zero;//zeros<dmat>(Ns,Ns);
}
  dmat Momentum::drhoim1(){
    // Todo: add this term!
    return zeros<dmat>(Ns,Ns);
}
  dmat Momentum::dUim1(){
    // Todo: add this term!
    return zeros<dmat>(Ns,Ns);
}

  dmat Momentum::dpip1(){
    TRACE(0,"Momentum::dpip1()");
    return tube.geom.Sf(i+1)*eye<dmat>(Ns,Ns)/(dxm+dxm);
  }
  dmat Momentum::dpim1(){
    return tube.geom.Sf(i-1)*eye<dmat>(Ns,Ns)/(dxm+dxm);
}
  Momentum::~Momentum(){
    TRACE(-5,"Momentum destructor");
}
} // namespace tube
