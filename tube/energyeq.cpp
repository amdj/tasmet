#include "energyeq.h"
#include "tubevertex.h"


#define ENERGY_SCALE  (vertex.gc->p0) //(1/v.gc.omg)
#define ENERGY_SCALE0 (1.0)//(1.0/pow(v.gc.M,2)) //(1/v.gc.omg)

namespace tube{
  Isentropic::Isentropic(TubeVertex& gp):TubeEquation(gp){
  }
  Isentropic::~Isentropic(){}
  vd Isentropic::Error(){
    TRACE(6,"Isentropic::Error()");
    TRACE(6,"Isentropic::Error() done");
    vd err(v.gc->Ns,fillwith::zeros);
    d T0=v.gc->T0;
    d p0=v.gc->p0;
    d rho0=v.gc->gas.rho(T0,p0);
    d gamma=v.gc->gas.gamma(T0);
    err+=(getp0()+v.p())/p0;
    err+=-1.0*v.gc->fDFT*pow(v.rho.tdata()/rho0,gamma);
    TRACE(6,"Isentropic::Error() done");
    return ENERGY_SCALE*err;
  }
  dmat Isentropic::dpi(){
    TRACE(1,"Isentropic::dpi()");
    dmat dpi(v.gc->Ns,v.gc->Ns,fillwith::eye);
    d p0=v.gc->p0;    
    dpi=dpi/p0;
    return ENERGY_SCALE*dpi;
  }
  dmat Isentropic::drhoi(){
    TRACE(1,"Isentropic::drhoi()"); 
    dmat drhoi(v.gc->Ns,v.gc->Ns,fillwith::zeros);

    d T0=v.gc->T0;
    d p0=v.gc->p0;
    d rho0=v.gc->gas.rho(T0,p0);

    d gamma=v.gc->gas.gamma(T0);
    drhoi+=-1.0*ENERGY_SCALE*(gamma/rho0)*v.gc->fDFT*
      diagmat(pow(v.rho.tdata()/rho0,(gamma-1.0)))*v.gc->iDFT;
    return drhoi;
  }
  Energy::Energy(TubeVertex& gp):TubeEquation(gp){
    TRACE(6,"Energy::Energy()"); 
    // Standard boundary condition is adiabatic-no-slip-wall
    if(v.i==0){			// Leftmost vertex
      TRACE(-1,"Leftmost vertex");

      
    } else if(v.i==v.Ncells-1){	// Rightmost vertex
      TRACE(-1,"Rightmost vertex");


    } else{			// Normal interior vertex

    }
    TRACE(0,"Energy constructor done");
  }
  vd Energy::Error(){		// Error in momentum equation
    TRACE(6,"Energy::Error()");
    vd error(v.gc->Ns,fillwith::zeros);
    d gamma=this->gamma();

    vd Uti=v.U.tdata();
    vd pti=v.p.tdata()+getp0t();
    vd Tti=v.T.tdata();
    error+=Wddt*v.gc->DDTfd*v.p()/(gamma-1.0);
    error+=Wgi*gamma*v.gc->fDFT*(pti%Uti)/(gamma-1.0);
    error+=Wji*v.gc->fDFT*(pti%Uti);
    error+=v.gc->fDFT*(Wc2*kappaL()%Tti+Wc3*kappaR()%Tti);
    
    if(v.i>0){
      vd Utim1=v.left->U.tdata();
      vd ptim1=v.left->p.tdata()+getp0t();
      vd Tim1=v.left->T.tdata();
      error+=Wgim1*gamma*v.gc->fDFT*(ptim1%Utim1)/(gamma-1.0);
      error+=Wjim1*v.gc->fDFT*(ptim1%Uti);
      error+=Wc1*v.gc->fDFT*(kappaL()%Tim1);
    }
    if(v.i<v.Ncells-1){
      vd Utip1=v.right->U.tdata();
      vd ptip1=v.right->p.tdata()+getp0t();
      vd Tip1=v.right->T.tdata();
      error+=Wgip1*gamma*v.gc->fDFT*(ptip1%Utip1)/(gamma-1.0);
      error+=Wjip1*v.gc->fDFT*(ptip1%Uti);
      error+=Wc4*v.gc->fDFT*(kappaR()%Tip1);
    }

    // (Boundary source term)
    TRACE(-1,"vertex.esource called from Vertex object..."<<v.esource());
    
    error+=v.esource();
    error(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*error;
  }

  dmat Energy::dpi(){
    TRACE(0,"Energy::dpi()");
    d T0=v.T(0);
    d gamma=v.gc->gas.gamma(T0);
    dmat dpi=zero;
    dmat diagUt=diagtmat(v.U);
    dpi+=(Wddt/(gamma-1.0))*v.gc->DDTfd;
    dpi+=Wgi*(gamma/(gamma-1.0))*v.gc->fDFT*diagUt*v.gc->iDFT;
    dpi+=Wji*v.gc->fDFT*diagUt*v.gc->iDFT;    
    dpi.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dpi;
  }

  dmat Energy::dUi(){
    TRACE(0,"Energy::dUi()");
    dmat dUi=zero;			    // Initialize with zeros
    d T0=v.T(0);
    d gamma=v.gc->gas.gamma(T0);
    dmat diagpt=diagmat(getp0t()+v.p.tdata());

    dUi+=Wgi*(gamma/(gamma-1.0))*v.gc->fDFT*diagpt*v.gc->iDFT;
    dUi+=Wji*v.gc->fDFT*diagpt*v.gc->iDFT;
    if(v.i>0)
      dUi+=Wjim1*v.gc->fDFT*diagmat(v.left->p.tdata()+getp0t())*v.gc->iDFT;
    if(v.i<v.Ncells-1)
      dUi+=Wjip1*v.gc->fDFT*diagmat(v.right->p.tdata()+getp0t())*v.gc->iDFT;
    dUi.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dUi;
  }
  dmat Energy::dTi(){
    TRACE(0,"Energy::dTi()");
    dmat dTi=zero;
    dTi+=v.gc->fDFT*(Wc2*diagmat(kappaL())+Wc3*diagmat(kappaR()))*v.gc->iDFT;
    dTi.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dTi;
  }
  dmat Energy::dpip1(){
    TRACE(0,"Energy::dpip1()");
    vd Uti=v.U.tdata();
    dmat dpip1=zero;
    d gamma=this->gamma();

    dpip1+=Wjip1*v.gc->fDFT*diagmat(Uti)*v.gc->iDFT;
    if(v.i<v.Ncells-1){
      vd Utip1=v.right->U.tdata();
      dpip1+=Wgip1*(gamma/(gamma-1))*v.gc->fDFT*diagmat(Utip1)*v.gc->iDFT;
    }
    dpip1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dpip1;
  }
  dmat Energy::dUip1(){
    TRACE(0,"Energy::dUip1()");
    d gamma=this->gamma();
    dmat dUip1=zero;
    if(v.i<v.Ncells-1){
      dUip1+=Wgip1*(gamma/(gamma-1.0))*
	v.gc->fDFT*diagmat(getp0t()+v.right->p.tdata())*v.gc->iDFT;
    }
    dUip1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dUip1;
  }
  dmat Energy::dTip1(){
    TRACE(0,"Energy::dTip1()");

    dmat dTip1=zero;
    if(v.i<v.Ncells-1)    
      dTip1+=Wc4*v.gc->fDFT*diagmat(kappaR())*v.gc->iDFT;
    dTip1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dTip1;
  }
  dmat Energy::dpim1(){
    TRACE(0,"Energy::dpim1()");
    d gamma=this->gamma();
    vd Uti=v.U.tdata();
    dmat dpim1=zero;
    dpim1+=Wjim1*v.gc->fDFT*diagmat(Uti)*v.gc->iDFT;
    if(v.i>0){
      vd Utim1=v.left->U.tdata();
      dpim1+=Wgim1*(gamma/(gamma-1))*v.gc->fDFT*diagmat(Utim1)*v.gc->iDFT;
    }
    dpim1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dpim1;
  }
  dmat Energy::dUim1(){
    TRACE(0,"Energy::dUim1()");
    d gamma=this->gamma();

    dmat dUim1=zero;
    if(v.i>0){
      dUim1+=Wgim1*(gamma/(gamma-1.0))*
	v.gc->fDFT*diagmat(getp0t()+v.left->p.tdata())*v.gc->iDFT;
    }
    dUim1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dUim1;
    }
  dmat Energy::dTim1(){
    TRACE(0,"Energy::dTim1()");
    dmat dTim1=zero;
    if(v.i>0)    
      dTim1+=Wc1*v.gc->fDFT*diagmat(kappaL())*v.gc->iDFT;
    dTim1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dTim1;
  }

  // From here on, auxiliary functions
  vd Energy::kappaL(){
    TRACE(0,"Energy::kappaL()");
    
    vd kappaL(v.gc->Ns,fillwith::zeros);
    vd Tti=v.T.tdata();
    vd kappait=v.gc->gas.kappa(Tti);

    if(v.i==0){
      vd Ttip1=v.right->T.tdata();
      vd kappaitp1=v.gc->gas.kappa(Ttip1);
      kappaL=v.wL1*kappaitp1+v.wL0*kappait;
    }
    else{
      vd Ttim1=v.left->T.tdata();
      vd kappaitm1=v.gc->gas.kappa(Ttim1);
      kappaL=v.wLr*kappait+v.wLl*kappaitm1;
      }
    // kappaL.zeros();		// WARNING
    return kappaL;
  }
  vd Energy::kappaR(){		// Returns thermal conductivity time domain data
    TRACE(0,"Energy::kappaR()");
    
    vd kappaR(v.gc->Ns,fillwith::zeros);
    vd Tti=v.T.tdata();
    vd kappait=v.gc->gas.kappa(Tti);

    if(v.i==v.Ncells-1){
      vd Ttim1=v.left->T.tdata();
      vd kappaitm1=v.gc->gas.kappa(Ttim1);
      kappaR=v.wRNm2*kappaitm1+v.wRNm1*kappait;
    }
    else{
      vd Ttip1=v.right->T.tdata();
      vd kappaitp1=v.gc->gas.kappa(Ttip1);
      kappaR=v.wRl*kappait+v.wRr*kappaitp1;
    }
    // kappaR.zeros();		// WARNING!!
    return kappaR;
  }
  d Energy::gamma(){
    d T0=v.T(0);
    return v.gc->gas.gamma(T0);
  }
  Energy::~Energy(){
    TRACE(-5,"Energy destructor");
  }

} // namespace tube
