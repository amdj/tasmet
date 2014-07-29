#include "energyeq.h"
#include "tubevertex.h"

// #define EN_VISCOSITY
#define ENERGY_SCALE  (1.0) //(1/v.gc.omg)
#define ENERGY_SCALE0 (1.0)//(1.0/pow(v.gc.M,2)) //(1/v.gc.omg)
// #define CONDUCTION


namespace tube{

  Energy::Energy(TubeVertex& i):TubeEquation(i){
    TRACE(6,"Energy::Energy()"); 
  }
  vd Energy::Error(){		// Error in momentum equation
    TRACE(6,"Energy::Error()");
    assert(v.gc!=NULL);
    vd error(v.gc->Ns,fillwith::zeros);
    d gamma=this->gamma();

    vd Uti=v.U.tdata();
    vd pti=v.p.tdata()+getp0t();
    vd Tti=v.T.tdata();
    error+=Wddt*v.gc->DDTfd*v.p()/(gamma-1.0);
    error+=Wgi*gamma*v.gc->fDFT*(pti%Uti)/(gamma-1.0);
    error+=Wji*v.gc->fDFT*(pti%Uti);
    #ifdef CONDUCTION
    error+=v.gc->fDFT*(Wc2*kappaL()%Tti+Wc3*kappaR()%Tti);
    #endif
    
    if(v.left!=NULL){
      vd Utim1=v.left->U.tdata();
      vd ptim1=v.left->p.tdata()+getp0t();
      vd Tim1=v.left->T.tdata();
      error+=Wgim1*gamma*v.gc->fDFT*(ptim1%Utim1)/(gamma-1.0);
      error+=Wjim1*v.gc->fDFT*(ptim1%Uti);
      #ifdef CONDUCTION
      error+=Wc1*v.gc->fDFT*(kappaL()%Tim1);
      #endif
    }
    if(v.right!=NULL){
      vd Utip1=v.right->U.tdata();
      vd ptip1=v.right->p.tdata()+getp0t();
      vd Tip1=v.right->T.tdata();
      error+=Wgip1*gamma*v.gc->fDFT*(ptip1%Utip1)/(gamma-1.0);
      error+=Wjip1*v.gc->fDFT*(ptip1%Uti);
      #ifdef CONDUCTION
      error+=Wc4*v.gc->fDFT*(kappaR()%Tip1);
      #endif
    }

    // Artificial viscosity terms      
    #ifdef EN_VISCOSITY
    const d& vVf=v.vVf;
    if(v.i>0 && v.i<v.Ncells-1){
      error+=-d_r()*(v.right->p() -v.p())*vVf;
      error+= d_l()*(v.p() -v.left->p())  *vVf;
    }
    else if(v.i==0){		// First v
      error+=-d_r()*(v.right->right->p()-v.right->p())*vVf;
      error+=d_l()*(v.right->p()-v.p())*vVf;
    }
    else {			// Last v
      error+=-d_r()*(v.p()-v.left->p())*vVf;
      error+=d_l()*(v.left->p()-v.left->left->p())*vVf;
    }
    #endif

    // (Boundary source term)
    // TRACE(-1,"vertex.esource called from Vertex object..."<<v.esource());
    
    error+=v.esource();
    return error;
  }

  dmat Energy::dpi(){
    TRACE(0,"Energy::dpi()");
    d T0=v.T(0);
    d gamma=v.gc->gas.gamma(T0);
    dmat dpi=zero;
    dmat diagUt=diagtmat(v.U);
    // Wddt = vVf, as defined in tubevertex.cpp
    dpi+=(Wddt/(gamma-1.0))*v.gc->DDTfd;
    dpi+=Wgi*(gamma/(gamma-1.0))*v.gc->fDFT*diagUt*v.gc->iDFT;
    dpi+=Wji*v.gc->fDFT*diagUt*v.gc->iDFT;    
    return dpi;
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
    // dUi.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dUi;
  }
  dmat Energy::dTi(){
    TRACE(0,"Energy::dTi()");
    dmat dTi=zero;
    #ifdef CONDUCTION
    dTi+=v.gc->fDFT*(Wc2*diagmat(kappaL())+Wc3*diagmat(kappaR()))*v.gc->iDFT;
    #endif
    // dTi.row(0)*=ENERGY_SCALE0;
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
    // dpip1.row(0)*=ENERGY_SCALE0;
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
    // dUip1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dUip1;
  }
  dmat Energy::dTip1(){
    TRACE(0,"Energy::dTip1()");

    dmat dTip1=zero;
    #ifdef CONDUCTION
    if(v.i<v.Ncells-1)    
      dTip1+=Wc4*v.gc->fDFT*diagmat(kappaR())*v.gc->iDFT;
    #endif
    // dTip1.row(0)*=ENERGY_SCALE0;
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
    // dpim1.row(0)*=ENERGY_SCALE0;
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
    // dUim1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dUim1;
  }
  dmat Energy::dTim1(){
    TRACE(0,"Energy::dTim1()");
    dmat dTim1=zero;
    #ifdef CONDUCTION
    if(v.i>0)    
      dTim1+=Wc1*v.gc->fDFT*diagmat(kappaL())*v.gc->iDFT;
    #endif
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
