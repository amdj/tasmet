#include "energyeq.h"
#include "vertex.h"
#include "tube.h"

#define ENERGY_SCALE (1.0) //(gc.p0) //(1/gc.omg)
#define ENERGY_SCALE0 (1.0)//(1.0/pow(gc.M,2)) //(1/gc.omg)

namespace tube{

  Energy::Energy(const Tube& tube,TubeVertex& gp):Equation(tube,gp){
    // Standard boundary condition is adiabatic-no-slip-wall
    if(i==0){			// Leftmost vertex
      TRACE(-1,"Leftmost vertex");

      Wgim1=0;
      Wgi=wRl;
      Wgip1=wRr;

      Wjim1=0;
      Wji=wL0-wRl;
      Wjip1=wL1-wRr;
      
      Wc1=0;
      Wc2=0;
      Wc3=SfR/dxp;
      Wc4=-SfR/dxp;
      
    } else if(i==Ncells-1){	// Rightmost vertex
      TRACE(-1,"Rightmost vertex");

      Wgim1=-wLl;
      Wgi=-wLr;
      Wgip1=0;

      Wjim1=wLl-wRNm2;
      Wji=wLr-wRNm1;
      Wjip1=0;

      Wc1=-SfL/dxm;
      Wc2=SfL/dxm;
      Wc3=0;
      Wc4=0;

    } else{			// Normal interior vertex

      Wgim1=-wLl;
      Wgi=wRl-wLr;
      Wgip1=wRr;

      Wjim1=wLl;
      Wji=wLr-wRl;
      Wjip1=-wRr;
      
      Wc1=-SfL/dxm;
      Wc2=SfL/dxm;
      Wc3=SfR/dxp;
      Wc4=-SfR/dxp;
    }
    TRACE(0,"Energy constructor done");
  }
  vd Energy::Error(){		// Error in momentum equation
    TRACE(0,"Energy::Error()");
    vd error(Ns,fillwith::zeros);
    d gamma=this->gamma();

    vd Uti=vertex.U.tdata();
    vd pti=vertex.p.tdata()+getp0t();
    vd Tti=vertex.T.tdata();
    error+=vVf*DDTfd*vertex.p()/(gamma-1.0);
    error+=Wgi*gamma*fDFT*(pti%Uti)/(gamma-1.0);
    error+=Wji*fDFT*(pti%Uti);
    error+=fDFT*(Wc2*kappaL()%Tti+Wc3*kappaR()%Tti);
    
    if(i>0){
      vd Utim1=left->U.tdata();
      vd ptim1=left->p.tdata()+getp0t();
      vd Tim1=left->T.tdata();
      error+=Wgim1*gamma*fDFT*(ptim1%Utim1)/(gamma-1.0);
      error+=Wjim1*fDFT*(ptim1%Uti);
      error+=Wc1*fDFT*(kappaL()%Tim1);
    }
    if(i<Ncells-1){
      vd Utip1=right->U.tdata();
      vd ptip1=right->p.tdata()+getp0t();
      vd Tip1=right->T.tdata();
      error+=Wgip1*gamma*fDFT*(ptip1%Utip1)/(gamma-1.0);
      error+=Wjip1*fDFT*(ptip1%Uti);
      error+=Wc4*fDFT*(kappaR()%Tip1);
    }

    // (Boundary source term)
    TRACE(-1,"vertex.esource called from Vertex object..."<<vertex.esource());
    
    error+=vertex.esource();
    error(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*error;
  }

  dmat Energy::dpi(){
    TRACE(0,"Energy::dpi()");
    d T0=vertex.T(0);
    d gamma=tube.gas.gamma(T0);
    dmat dpi=zero;
    dmat diagUt=diagtmat(vertex.U);
    dpi+=(vVf/(gamma-1.0))*DDTfd;
    dpi+=Wgi*(gamma/(gamma-1.0))*fDFT*diagUt*iDFT;
    dpi+=Wji*fDFT*diagUt*iDFT;    
    dpi.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dpi;
  }

  dmat Energy::dUi(){
    TRACE(0,"Energy::dUi()");
    dmat dUi=zero;			    // Initialize with zeros
    d T0=vertex.T(0);
    d gamma=tube.gas.gamma(T0);
    dmat diagpt=diagmat(getp0t()+vertex.p.tdata());

    dUi+=Wgi*(gamma/(gamma-1.0))*fDFT*diagpt*iDFT;
    dUi+=Wji*fDFT*diagpt*iDFT;
    if(i>0)
      dUi+=Wjim1*fDFT*diagmat(left->p.tdata()+getp0t())*iDFT;
    if(i<Ncells-1)
      dUi+=Wjip1*fDFT*diagmat(right->p.tdata()+getp0t())*iDFT;
    dUi.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dUi;
  }
  dmat Energy::dTi(){
    TRACE(0,"Energy::dTi()");
    dmat dTi=zero;
    dTi+=fDFT*(Wc2*diagmat(kappaL())+Wc3*diagmat(kappaR()))*iDFT;
    dTi.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dTi;
  }
  dmat Energy::dpip1(){
    TRACE(0,"Energy::dpip1()");
    vd Uti=vertex.U.tdata();
    dmat dpip1=zero;
    d gamma=this->gamma();

    dpip1+=Wjip1*fDFT*diagmat(Uti)*iDFT;
    if(i<Ncells-1){
      vd Utip1=right->U.tdata();
      dpip1+=Wgip1*(gamma/(gamma-1))*fDFT*diagmat(Utip1)*iDFT;
    }
    dpip1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dpip1;
  }
  dmat Energy::dUip1(){
    TRACE(0,"Energy::dUip1()");
    d gamma=this->gamma();
    dmat dUip1=zero;
    if(i<Ncells-1){
      dUip1+=Wgip1*(gamma/(gamma-1.0))*fDFT*diagmat(getp0t()+right->p.tdata())*iDFT;
    }
    dUip1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dUip1;
  }
  dmat Energy::dTip1(){
    TRACE(0,"Energy::dTip1()");
    dmat dTip1=zero;
    if(i<Ncells-1)    
      dTip1+=Wc4*fDFT*diagmat(kappaR())*iDFT;
    dTip1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dTip1;
  }
  dmat Energy::dpim1(){
    TRACE(0,"Energy::dpim1()");
    d gamma=this->gamma();

    vd Uti=vertex.U.tdata();
    dmat dpim1=zero;
    dpim1+=Wjim1*fDFT*diagmat(Uti)*iDFT;
    if(i>0){
      vd Utim1=left->U.tdata();
      dpim1+=Wgim1*(gamma/(gamma-1))*fDFT*diagmat(Utim1)*iDFT;
    }
    dpim1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dpim1;
  }
  dmat Energy::dUim1(){
    TRACE(0,"Energy::dUim1()");
    d gamma=this->gamma();

    dmat dUim1=zero;
    if(i>0){
      dUim1+=Wgim1*(gamma/(gamma-1.0))*fDFT*diagmat(getp0t()+left->p.tdata())*iDFT;
    }
    dUim1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dUim1;
    }
  dmat Energy::dTim1(){
    TRACE(0,"Energy::dTim1()");
    dmat dTim1=zero;
    if(i>0)    
      dTim1+=Wc1*fDFT*diagmat(kappaL())*iDFT;
    dTim1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dTim1;
  }

  // From here on, auxiliary functions
  vd Energy::kappaL(){
    TRACE(0,"Energy::kappaL()");
    
    vd kappaL(Ns,fillwith::zeros);
    vd Tti=vertex.T.tdata();
    vd kappait=tube.gas.kappa(Tti);

    if(i==0){
      vd Ttip1=right->T.tdata();
      vd kappaitp1=tube.gas.kappa(Ttip1);
      kappaL=wL1*kappaitp1+wL0*kappait;
    }
    else{
      vd Ttim1=left->T.tdata();
      vd kappaitm1=tube.gas.kappa(Ttim1);
      kappaL=wLr*kappait+wLl*kappaitm1;
      }
    // kappaL.zeros();		// WARNING
    return kappaL;
  }
  vd Energy::kappaR(){		// Returns thermal conductivity time domain data
    TRACE(0,"Energy::kappaR()");
    
    vd kappaR(Ns,fillwith::zeros);
    vd Tti=vertex.T.tdata();
    vd kappait=tube.gas.kappa(Tti);

    if(i==Ncells-1){
      vd Ttim1=left->T.tdata();
      vd kappaitm1=tube.gas.kappa(Ttim1);
      kappaR=wRNm2*kappaitm1+wRNm1*kappait;
    }
    else{
      vd Ttip1=right->T.tdata();
      vd kappaitp1=tube.gas.kappa(Ttip1);
      kappaR=wRl*kappait+wRr*kappaitp1;
    }
    // kappaR.zeros();		// WARNING!!
    return kappaR;
  }
  d Energy::gamma(){
    d T0=vertex.T(0);
    return tube.gas.gamma(T0);
  }
  Energy::~Energy(){
    TRACE(-5,"Energy destructor");
  }

  Isentropic::Isentropic(const Tube& tube,TubeVertex& gp):Equation(tube,gp){
  }
  Isentropic::~Isentropic(){}
  vd Isentropic::Error(){
    TRACE(0,"Isentropic::Error()");
    vd err(Ns,fillwith::zeros);
    d T0=tube.gc.T0;
    d p0=tube.gc.p0;
    d gamma=tube.gas.gamma(T0);
    err+=(getp0()+vertex.p())/p0;
    err+=-1.0*fDFT*pow(vertex.T.tdata()/T0,gamma/(gamma-1.0));
    return ENERGY_SCALE*err;
  }
  dmat Isentropic::dpi(){
    TRACE(0,"Isentropic::dpi()");
    dmat dpi(Ns,Ns,fillwith::eye);
    d p0=tube.gc.p0;    
    dpi=dpi/p0;
    return ENERGY_SCALE*dpi;
  }
  dmat Isentropic::dTi(){
    TRACE(0,"Isentropic::dTi()");    
    dmat dTi(Ns,Ns,fillwith::zeros);
    d T0=tube.gc.T0;    
    d gamma=tube.gas.gamma(T0);
    dTi+=-1.0*ENERGY_SCALE*(gamma/((gamma-1.0)*T0))*fDFT*diagmat(pow(vertex.T.tdata()/T0,-1.0/gamma))*iDFT;
    return dTi;
  }
} // namespace tube




