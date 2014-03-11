#include "continuityeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{


  Continuity::Continuity(Tube* tube,TubeVertex* tgp):Equation(tube,tgp){
    TRACE(0,"Continuity constructor done");
    TRACE(-2,"rho:"<< vertex->rho());
  }
  dmat Continuity::operator()(){
    TRACE(0,"Continuity::operator()");
    return Equation::operator()();
  }
  vd Continuity::Error(){
    TRACE(0,"Continuity::Error()");
    TRACE(-1,"i: "<<i);
    assert(i>0 && i<tube->geom.gp-1);
    vd error(Ns,fillwith::zeros);
    error+=Sf*vop.DDTfd*vertex->rho();
    error=error+vop.fDFT*(tube->gps[i+1]->U.tdata()%tube->gps[i+1]->rho.tdata())/(dxp+dxm);
    error=error-1.0*vop.fDFT*(tube->gps[i-1]->U.tdata()%tube->gps[i-1]->rho.tdata())/(dxp+dxm);
    return error;
  }
  dmat  Continuity::drhoim1(){
    dmat Uim1td(Ns,Ns,fillwith::zeros);			// Ui-1 in time domain
    Uim1td.diag()=tube->gps.at(i-1)->U.tdata();
    dmat rhoim1=-1.0/(dxp+dxm)*vop.fDFT*Uim1td*vop.iDFT;
    return rhoim1;
  }
  dmat Continuity::dUim1(){
    TRACE(0,"Continuity::dUim1()");
    dmat rhoim1td(Ns,Ns,fillwith::zeros);// rhoi-1 in time domain
    rhoim1td.diag()=tube->gps.at(i-1)->rho.tdata();
    dmat dUim1=-1.0/(dxp+dxm)*vop.fDFT*rhoim1td*vop.iDFT;
    return dUim1;
  }
  dmat Continuity::drhoip1(){
    dmat Uip1td(Ns,Ns,fillwith::zeros);	      	// Ui+1 in time domain
    Uip1td.diag()=tube->gps.at(i+1)->U.tdata();
    dmat drhoip1=1.0/(dxp+dxm)*vop.fDFT*Uip1td*vop.iDFT;
    return drhoip1;
  }
  dmat Continuity::dUip1(){
    dmat rhoip1td(Ns,Ns,fillwith::zeros); 	// rhoi+1 in time domain
    rhoip1td.diag()=tube->gps.at(i+1)->rho.tdata();
    dmat dUip1=1.0/(dxp+dxm)*vop.fDFT*rhoip1td*vop.iDFT;
    return dUip1;
  }
  dmat Continuity::drhoi(){
    TRACE(0,"Continuity::drhoi()");
    return vop.DDTfd;
  }
  Continuity::~Continuity(){}

} // Namespace tube
