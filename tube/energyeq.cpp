#include "energyeq.h"
#include "tubevertex.h"


#define ENERGY_SCALE  (gc->p0) //(1/gc.omg)
#define ENERGY_SCALE0 (1.0)//(1.0/pow(gc.M,2)) //(1/gc.omg)

namespace tube{
  Isentropic::Isentropic(TubeVertex& gp):TubeEquation(gp){
  }
  Isentropic::~Isentropic(){}
  vd Isentropic::Error(){
    TRACE(6,"Isentropic::Error()");
    vd err(gc->Ns,fillwith::zeros);
    const Globalconf& gc1=*this->gc;
    d T0=gc1.T0;
    d p0=gc1.p0;
    d rho0=gc1.gas.rho(T0,p0);
    d gamma=gc1.gas.gamma(T0);
    err+=(getp0()+vertex.p())/p0;
    err+=-1.0*gc->fDFT*pow(vertex.rho.tdata()/rho0,gamma);
    return ENERGY_SCALE*err;
  }
  dmat Isentropic::dpi(){
    TRACE(1,"Isentropic::dpi()");
    dmat dpi(gc->Ns,gc->Ns,fillwith::eye);
    d p0=gc->p0;    
    dpi=dpi/p0;
    return ENERGY_SCALE*dpi;
  }
  dmat Isentropic::drhoi(){
    TRACE(1,"Isentropic::drhoi()"); 
    dmat drhoi(gc->Ns,gc->Ns,fillwith::zeros);

    d T0=gc->T0;
    d p0=gc->p0;
    d rho0=gc->gas.rho(T0,p0);

    d gamma=gc->gas.gamma(T0);
    drhoi+=-1.0*ENERGY_SCALE*(gamma/rho0)*gc->fDFT*
      diagmat(pow(vertex.rho.tdata()/rho0,(gamma-1.0)))*gc->iDFT;
    return drhoi;
  }
  Energy::Energy(TubeVertex& gp):TubeEquation(gp){
    TRACE(6,"Energy::Energy()"); 
    // Standard boundary condition is adiabatic-no-slip-wall
    if(i==0){			// Leftmost vertex
      TRACE(-1,"Leftmost vertex");

      
    } else if(i==Ncells-1){	// Rightmost vertex
      TRACE(-1,"Rightmost vertex");


    } else{			// Normal interior vertex

    }
    TRACE(0,"Energy constructor done");
  }
  vd Energy::Error(){		// Error in momentum equation
    TRACE(6,"Energy::Error()");
    vd error(gc->Ns,fillwith::zeros);
    d gamma=this->gamma();

    vd Uti=vertex.U.tdata();
    vd pti=vertex.p.tdata()+getp0t();
    vd Tti=vertex.T.tdata();
    error+=Wddt*gc->DDTfd*vertex.p()/(gamma-1.0);
    error+=Wgi*gamma*gc->fDFT*(pti%Uti)/(gamma-1.0);
    error+=Wji*gc->fDFT*(pti%Uti);
    error+=gc->fDFT*(Wc2*kappaL()%Tti+Wc3*kappaR()%Tti);
    
    if(i>0){
      vd Utim1=left->U.tdata();
      vd ptim1=left->p.tdata()+getp0t();
      vd Tim1=left->T.tdata();
      error+=Wgim1*gamma*gc->fDFT*(ptim1%Utim1)/(gamma-1.0);
      error+=Wjim1*gc->fDFT*(ptim1%Uti);
      error+=Wc1*gc->fDFT*(kappaL()%Tim1);
    }
    if(i<Ncells-1){
      vd Utip1=right->U.tdata();
      vd ptip1=right->p.tdata()+getp0t();
      vd Tip1=right->T.tdata();
      error+=Wgip1*gamma*gc->fDFT*(ptip1%Utip1)/(gamma-1.0);
      error+=Wjip1*gc->fDFT*(ptip1%Uti);
      error+=Wc4*gc->fDFT*(kappaR()%Tip1);
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
    d gamma=gc->gas.gamma(T0);
    dmat dpi=zero;
    dmat diagUt=diagtmat(vertex.U);
    dpi+=(Wddt/(gamma-1.0))*gc->DDTfd;
    dpi+=Wgi*(gamma/(gamma-1.0))*gc->fDFT*diagUt*gc->iDFT;
    dpi+=Wji*gc->fDFT*diagUt*gc->iDFT;    
    dpi.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dpi;
  }

  dmat Energy::dUi(){
    TRACE(0,"Energy::dUi()");
    dmat dUi=zero;			    // Initialize with zeros
    d T0=vertex.T(0);
    d gamma=gc->gas.gamma(T0);
    dmat diagpt=diagmat(getp0t()+vertex.p.tdata());

    dUi+=Wgi*(gamma/(gamma-1.0))*gc->fDFT*diagpt*gc->iDFT;
    dUi+=Wji*gc->fDFT*diagpt*gc->iDFT;
    if(i>0)
      dUi+=Wjim1*gc->fDFT*diagmat(left->p.tdata()+getp0t())*gc->iDFT;
    if(i<Ncells-1)
      dUi+=Wjip1*gc->fDFT*diagmat(right->p.tdata()+getp0t())*gc->iDFT;
    dUi.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dUi;
  }
  dmat Energy::dTi(){
    TRACE(0,"Energy::dTi()");
    dmat dTi=zero;
    dTi+=gc->fDFT*(Wc2*diagmat(kappaL())+Wc3*diagmat(kappaR()))*gc->iDFT;
    dTi.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dTi;
  }
  dmat Energy::dpip1(){
    TRACE(0,"Energy::dpip1()");
    vd Uti=vertex.U.tdata();
    dmat dpip1=zero;
    d gamma=this->gamma();

    dpip1+=Wjip1*gc->fDFT*diagmat(Uti)*gc->iDFT;
    if(i<Ncells-1){
      vd Utip1=right->U.tdata();
      dpip1+=Wgip1*(gamma/(gamma-1))*gc->fDFT*diagmat(Utip1)*gc->iDFT;
    }
    dpip1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dpip1;
  }
  dmat Energy::dUip1(){
    TRACE(0,"Energy::dUip1()");
    d gamma=this->gamma();
    dmat dUip1=zero;
    if(i<Ncells-1){
      dUip1+=Wgip1*(gamma/(gamma-1.0))*
	gc->fDFT*diagmat(getp0t()+right->p.tdata())*gc->iDFT;
    }
    dUip1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dUip1;
  }
  dmat Energy::dTip1(){
    TRACE(0,"Energy::dTip1()");
    dmat dTip1=zero;
    if(i<Ncells-1)    
      dTip1+=Wc4*gc->fDFT*diagmat(kappaR())*gc->iDFT;
    dTip1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dTip1;
  }
  dmat Energy::dpim1(){
    TRACE(0,"Energy::dpim1()");
    d gamma=this->gamma();

    vd Uti=vertex.U.tdata();
    dmat dpim1=zero;
    dpim1+=Wjim1*gc->fDFT*diagmat(Uti)*gc->iDFT;
    if(i>0){
      vd Utim1=left->U.tdata();
      dpim1+=Wgim1*(gamma/(gamma-1))*gc->fDFT*diagmat(Utim1)*gc->iDFT;
    }
    dpim1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dpim1;
  }
  dmat Energy::dUim1(){
    TRACE(0,"Energy::dUim1()");
    d gamma=this->gamma();

    dmat dUim1=zero;
    if(i>0){
      dUim1+=Wgim1*(gamma/(gamma-1.0))*
	gc->fDFT*diagmat(getp0t()+left->p.tdata())*gc->iDFT;
    }
    dUim1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dUim1;
    }
  dmat Energy::dTim1(){
    TRACE(0,"Energy::dTim1()");
    dmat dTim1=zero;
    if(i>0)    
      dTim1+=Wc1*gc->fDFT*diagmat(kappaL())*gc->iDFT;
    dTim1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dTim1;
  }

  // From here on, auxiliary functions
  vd Energy::kappaL(){
    TRACE(0,"Energy::kappaL()");
    
    vd kappaL(gc->Ns,fillwith::zeros);
    vd Tti=vertex.T.tdata();
    vd kappait=gc->gas.kappa(Tti);

    if(i==0){
      vd Ttip1=right->T.tdata();
      vd kappaitp1=gc->gas.kappa(Ttip1);
      kappaL=vertex.wL1*kappaitp1+vertex.wL0*kappait;
    }
    else{
      vd Ttim1=left->T.tdata();
      vd kappaitm1=gc->gas.kappa(Ttim1);
      kappaL=vertex.wLr*kappait+vertex.wLl*kappaitm1;
      }
    // kappaL.zeros();		// WARNING
    return kappaL;
  }
  vd Energy::kappaR(){		// Returns thermal conductivity time domain data
    TRACE(0,"Energy::kappaR()");
    
    vd kappaR(gc->Ns,fillwith::zeros);
    vd Tti=vertex.T.tdata();
    vd kappait=gc->gas.kappa(Tti);

    if(i==Ncells-1){
      vd Ttim1=left->T.tdata();
      vd kappaitm1=gc->gas.kappa(Ttim1);
      kappaR=vertex.wRNm2*kappaitm1+vertex.wRNm1*kappait;
    }
    else{
      vd Ttip1=right->T.tdata();
      vd kappaitp1=gc->gas.kappa(Ttip1);
      kappaR=vertex.wRl*kappait+vertex.wRr*kappaitp1;
    }
    // kappaR.zeros();		// WARNING!!
    return kappaR;
  }
  d Energy::gamma(){
    d T0=vertex.T(0);
    return gc->gas.gamma(T0);
  }
  Energy::~Energy(){
    TRACE(-5,"Energy destructor");
  }

} // namespace tube
