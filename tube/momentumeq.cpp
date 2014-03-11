#include "momentumeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{

  Momentum::Momentum(Tube* tube,TubeVertex* gp):Equation(tube,gp){
    TRACE(0,"Momentum constructor done");
  }
  dmat Momentum::operator()(){
    TRACE(0,"Momentum::operator()");
    return Equation::operator()();
  }
  vd Momentum::Error(){		// Error in momentum equation
    TRACE(0,"Momentum::Error()");
    vd error(Ns,fillwith::zeros);
    error+=vop.DDTfd*vop.fDFT*(vertex->U.tdata()%vertex->rho.tdata());
    error+=tube->gps.at(i+1)->p()*tube->geom.Sf(i+1)/(dxp+dxm);
    error+=-1.0*tube->gps.at(i-1)->p()*tube->geom.Sf(i-1)/(dxp+dxm);
    error+=tube->drag(i);
    vd rhotip1=tube->gps[i+1]->rho.tdata();
    vd rhotim1=tube->gps[i-1]->rho.tdata();
    vd Uip1=tube->gps[i+1]->U.tdata();
    vd Uim1=tube->gps[i-1]->U.tdata();
    error+=vop.fDFT*(rhotip1%Uip1%Uip1)/(tube->geom.Sf(i+1)*(dxp+dxm));
    error+=-1.0*vop.fDFT*(rhotim1%Uim1%Uim1)/(tube->geom.Sf(i-1)*(dxp+dxm));
    return error;
  }
  dmat Momentum::drhoi(){
    TRACE(0,"Momentum::drhoi()");
    dmat Utdiag(Ns,Ns,fillwith::zeros);
    Utdiag.diag()=tube->gps.at(i)->U.tdata();
    return vop.DDTfd*vop.fDFT*Utdiag*vop.iDFT;
  }
  dmat Momentum::dUi(){
    TRACE(0,"Momentum::dUi()");
    TRACE(-1,"Ns:"<<Ns);
    dmat rhotdiag(Ns,Ns,fillwith::zeros);
    rhotdiag.diag()=tube->gps.at(i)->rho.tdata();
    TRACE(-1,"rhotdiag:"<<rhotdiag);
    dmat ddragdU=tube->drag.dUi(i);
    return ddragdU+vop.DDTfd*vop.fDFT*rhotdiag*vop.iDFT;
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
    return tube->geom.Sf(i+1)*eye<dmat>(Ns,Ns)/(dxm+dxm);
  }
  dmat Momentum::dpim1(){
    return tube->geom.Sf(i-1)*eye<dmat>(Ns,Ns)/(dxm+dxm);
}
  Momentum::~Momentum(){
    TRACE(-5,"Momentum destructor");
}
} // namespace tube
