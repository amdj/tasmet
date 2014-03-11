#include "energyeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{

  Energy::Energy(Tube* tube,TubeVertex* gp):Equation(tube,gp){
    TRACE(0,"Energy constructor done");
  }
  vd Energy::Error(){		// Error in momentum equation
    TRACE(0,"Energy::Error()");
    vd error(Ns,fillwith::zeros);
    d T0=vertex->T(0);
    d gamma=tube->gas.gamma(T0);
    TRACE(-1,"gamma: " << gamma);
    error+=tube->geom.Sf(i)*vop.DDTfd*vertex->p()/(gamma-1.0); // dpdt term

    vd ptip1=tube->gps[i+1]->p.tdata();
    vd ptim1=tube->gps[i-1]->p.tdata();

    vd Utip1=tube->gps[i+1]->U.tdata();
    vd Utim1=tube->gps[i-1]->U.tdata();
    vd Uti=vertex->U.tdata();
    error+=(gamma/(gamma-1.0))*vop.fDFT*(Utip1%ptip1)/(dxp+dxm); // Second term
    error+=-1.0*(gamma/(gamma-1.0))*vop.fDFT*(Utim1%ptim1)/(dxp+dxm); // Third term
    error+=-1.0*vop.fDFT*(Uti%(ptip1/(dxp+dxm)-ptim1/(dxp+dxm))); // Last term left hand side

    return error;
  }
  dmat Energy::dpi(){
    TRACE(0,"Energy::dpi()");
    d T0=vertex->T(0);
    d gamma=tube->gas.gamma(T0);
    return vertex->vSf*vop.DDTfd/(gamma-1.0);
  }
  dmat Energy::operator()(){
    TRACE(0,"Energy::operator()");
    return Equation::operator()();
  }
  dmat Energy::dUi(){
    TRACE(0,"Energy::dUi()");
    // dmat Utdiag(Ns,Ns,fillwith::zeros);
    // dmat result(Ns,Ns,fillwith::zeros);
    dmat dpdxtddiag(Ns,Ns,fillwith::zeros); // The derivative of the pressure in time domain
    vd ptip1=tube->gps[i+1]->p.tdata();
    vd ptim1=tube->gps[i-1]->p.tdata();
    dpdxtddiag.diag()=(ptip1-ptim1)/(dxp+dxm);
    return -1.0*vop.fDFT*dpdxtddiag*vop.iDFT;
  }
  dmat Energy::dUip1(){
    TRACE(0,"Energy::dUip1()");
    d T0=vertex->T(0);
    d gamma=tube->gas.gamma(T0);
    dmat pip1tdiag(Ns,Ns,fillwith::zeros); // The derivative of the pressure in time domain
    vd ptip1=tube->gps[i+1]->p.tdata();
    pip1tdiag.diag()=ptip1;
    return ((gamma-1.0)/(gamma*(dxp+dxm)))*vop.fDFT*pip1tdiag*vop.iDFT;
  }
  dmat Energy::dUim1(){
    TRACE(0,"Energy::dUim1()");
    d T0=vertex->T(0);
    d gamma=tube->gas.gamma(T0);
    dmat pim1tdiag(Ns,Ns,fillwith::zeros); // The derivative of the pressure in time domain
    vd ptim1=tube->gps[i+1]->p.tdata();
    pim1tdiag.diag()=ptim1;
    return ((gamma-1.0)/(gamma*(dxp+dxm)))*vop.fDFT*pim1tdiag*vop.iDFT;
  }
  dmat Energy::dpip1(){
    TRACE(0,"Energy::dpip1()");
    dmat term1(Ns,Ns,fillwith::zeros);
    dmat term2(Ns,Ns,fillwith::zeros);
    d T0=vertex->T(0);
    d gamma=tube->gas.gamma(T0);
    vd Utip1=tube->gps[i+1]->U.tdata();
    vd Uti=vertex->U.tdata();

    term1.diag()=(gamma-1.0)*Utip1/gamma;
    term2.diag()=Uti;
    return (1.0/(dxp+dxm))*vop.fDFT*(term1-term2)*vop.iDFT;
  }
  dmat Energy::dpim1(){
    TRACE(0,"Energy::dpim1()");
    dmat term1(Ns,Ns,fillwith::zeros);
    dmat term2(Ns,Ns,fillwith::zeros);
    d T0=vertex->T(0);
    d gamma=tube->gas.gamma(T0);
    vd Utim1=tube->gps[i-1]->U.tdata();
    vd Uti=vertex->U.tdata();

    term1.diag()=(gamma-1.0)*Utim1/gamma;
    term2.diag()=Uti;
    return (-1.0/(dxp+dxm))*vop.fDFT*(term1-term2)*vop.iDFT;
  }
  Energy::~Energy(){
    TRACE(-5,"Energy destructor");
  }
} // namespace tube
