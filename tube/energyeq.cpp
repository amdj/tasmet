#include "energyeq.h"
#include "tubevertex.h"
#include "conduction.h"

// #define ENERGY_SCALE (1.0/(v.gc->p0))
#define ENERGY_SCALE (1.0)
#define EN_VISCOSITY


namespace tube{

  Energy::Energy(TubeVertex& i):TubeEquation(i){
    TRACE(6,"Energy::Energy()"); 
  }
  void Energy::show(){
    cout << "--------- Showing energy weight factors for i=" << v.i <<"\n" \
	 << "Wddt      : "<<Wddt      <<"\n"					\
	 << "Wgim1     : "<<Wgim1<<"\n"						\
	 << "Wgi       : "<<Wgi      <<"\n"					\
	 << "Wgip1     : "<<Wgip1<<"\n"						\
	 << "Wjim1     : "<<Wjim1      <<"\n"					\
	 << "Wji       : "<<Wji      <<"\n"					\
	 << "Wjip1     : "<<Wjip1      <<"\n"					\
      ;      

  }
  vd Energy::Error(){		// Error in momentum equation
    TRACE(6,"Energy::Error()");
    assert(v.gc!=NULL);
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    vd error(v.gc->Ns,fillwith::zeros);
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    vd Uti=v.U.tdata();
    vd pti=v.p.tdata()+getp0t();
    vd Tti=v.T.tdata();
    error+=Wddt*DDTfd*v.p()/(gamma-1.0);
    error+=Wgi*gamfac*fDFT*(pti%Uti);
    error+=Wji*fDFT*(pti%Uti);
    error+=fDFT*((Wc2*kappaL()+Wc3*kappaR())%Tti);
    
    if(v.left!=NULL){
      vd Utim1=v.left->U.tdata();
      vd ptim1=v.left->p.tdata()+getp0t();
      vd Tim1=v.left->T.tdata();
      error+=Wgim1*gamfac*fDFT*(ptim1%Utim1);
      error+=Wjim1*fDFT*(ptim1%Uti);
      // Conduction term
      error+=Wc1*fDFT*(kappaL()%Tim1);

    }
    if(v.right!=NULL){
      vd Utip1=v.right->U.tdata();
      vd ptip1=v.right->p.tdata()+getp0t();
      vd Tip1=v.right->T.tdata();
      error+=Wgip1*gamfac*fDFT*(ptip1%Utip1);
      error+=Wjip1*fDFT*(ptip1%Uti);
      error+=Wc4*fDFT*(kappaR()%Tip1);
    }

    // Artificial viscosity terms      
    #ifdef EN_VISCOSITY
    const d& vSf=v.vSf;
    if(v.i>0 && v.i<v.Ncells-1){
      error+=-d_r()*(v.right->p() -v.p())*vSf;
      error+= d_l()*(v.p() -v.left->p())  *vSf;
    }
    else if(v.i==0){		// First v
      error+=-d_r()*(v.right->right->p()-v.right->p())*vSf;
      error+=d_l()*(v.right->p()-v.p())*vSf;
    }
    else {			// Last v
      error+=-d_r()*(v.p()-v.left->p())*vSf;
      error+=d_l()*(v.left->p()-v.left->left->p())*vSf;
    }
    #endif
    // (Boundary source term)
    error+=v.esource();
    return ENERGY_SCALE*error;
  }
  dmat Energy::dpi(){
    TRACE(0,"Energy::dpi()");
    d T0=v.T(0);
    d gamma=this->gamma();
    dmat dpi=zero;
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    // Wddt = vVf, as defined in tubevertex.cpp
    d gamfac=gamma/(gamma-1.0);
    dpi+=(Wddt/(gamma-1.0))*DDTfd;
    dpi+=Wgi*gamfac*fDFT*v.U.diagt()*iDFT;
    dpi+=Wji*fDFT*v.U.diagt()*iDFT;
    // Artificial viscosity terms    
    #ifdef EN_VISCOSITY
    const d& vVf=v.vVf;
    const d& vSf=v.vSf;
    if(v.i>0 && v.i<v.Ncells-1){
      dpi+=(d_l()+d_r())*vSf;	// Middle vertex
    }
    else if(v.i==0)
      dpi+=-d_l()*vSf;	// First vertex
    else		
      dpi+=-d_r()*vSf;	// Last vertex
    #endif
    return ENERGY_SCALE*dpi;
  }
  dmat Energy::dpim1(){
    TRACE(0,"Energy::dpim1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);    
    vd Uti=v.U.tdata();
    dmat dpim1=zero;
    if(v.left!=NULL){
      dpim1+=Wgim1*gamfac*fDFT*v.left->U.diagt()*iDFT;
    }
    dpim1+=Wjim1*fDFT*v.U.diagt()*iDFT;

    // Artificial viscosity terms
    #ifdef EN_VISCOSITY
    const d& vSf=v.vSf;
    if(v.i>0 && v.i<v.Ncells-1){
      dpim1+=-d_l()*vSf;
    }
    else if(v.i==v.Ncells-1){		// Last vertex
      dpim1+=(d_l()+d_r())*vSf;
    }
    #endif
    return ENERGY_SCALE*dpim1;
  }
  dmat Energy::dpip1(){
    TRACE(0,"Energy::dpip1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    vd Uti=v.U.tdata();
    dmat dpip1=zero;
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    if(v.right!=NULL){
      dpip1+=Wgip1*gamfac*fDFT*v.right->U.diagt()*iDFT;
    }
    dpip1+=Wjip1*fDFT*v.U.diagt()*iDFT;
    // Artificial viscosity terms
    #ifdef EN_VISCOSITY    
    const d& vSf=v.vSf;
    if(v.i>0 && v.i<v.Ncells-1){
      dpip1+=-d_r()*vSf;
    }
    else if(v.i==0){		// First vertex
      dpip1+=(d_l()+d_r())*vSf;
    }
    #endif
    return ENERGY_SCALE*dpip1;
  }
  dmat Energy::dUi(){
    TRACE(0,"Energy::dUi()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    dmat dUi=zero;			    // Initialize with zeros
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    dmat diagpt=diagmat(getp0t()+v.p.tdata());
    dUi+=Wgi*gamfac*fDFT*diagpt*iDFT;
    dUi+=Wji*fDFT*diagpt*iDFT;
    if(v.left!=NULL)
      dUi+=Wjim1*fDFT*diagmat(v.left->p.tdata()+getp0t())*iDFT;
    if(v.right!=NULL)
      dUi+=Wjip1*fDFT*diagmat(v.right->p.tdata()+getp0t())*iDFT;
    return ENERGY_SCALE*dUi;
  }
  dmat Energy::dUim1(){
    TRACE(0,"Energy::dUim1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    dmat dUim1=zero;
    if(v.left!=NULL){
      dUim1+=Wgim1*gamfac*fDFT*diagmat(getp0t()+v.left->p.tdata())*iDFT;
    }
    return ENERGY_SCALE*dUim1;
  }
  dmat Energy::dUip1(){
    TRACE(0,"Energy::dUip1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    dmat dUip1=zero;
    if(v.right!=NULL){
      dUip1+=Wgip1*gamfac*fDFT*diagmat(getp0t()+v.right->p.tdata())*iDFT;
    }
    return ENERGY_SCALE*dUip1;
  }
  dmat Energy::dpip2(){
    dmat dpip2=zero;
    #ifdef EN_VISCOSITY
    if(v.i==0 && v.left==NULL){
      const d& vSf=v.vSf;
      dpip2+=-d_r()*vSf;
    }
    #endif 
    return ENERGY_SCALE*dpip2;
  }
  dmat Energy::dpim2(){
    dmat dpim2=zero;
    #ifdef EN_VISCOSITY
    if(v.i==v.Ncells-1 && v.right==NULL){
      const d& vSf=v.vSf;
      dpim2+=-d_l()*vSf;
    }
    #endif 
    return ENERGY_SCALE*dpim2;
  }
  // ############################## TEMPERATURE AND CONDUCTION TERMS
  dmat Energy::dTip1(){
    TRACE(0,"Energy::dTip1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      

    dmat dTip1=zero;
    if(v.right!=NULL)    
      dTip1+=Wc4*fDFT*diagmat(kappaR())*iDFT;
    // dTip1.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dTip1;
  }
  dmat Energy::dTi(){
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    TRACE(0,"Energy::dTi()");
    dmat dTi=zero;
    dTi+=fDFT*diagmat(Wc2*kappaL()+Wc3*kappaR())*iDFT;
    // dTi.row(0)*=ENERGY_SCALE0;
    return ENERGY_SCALE*dTi;
  }
  dmat Energy::dTim1(){
    TRACE(0,"Energy::dTim1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    dmat dTim1=zero;
    if(v.left!=NULL)    
      dTim1+=Wc1*fDFT*diagmat(kappaL())*iDFT;
    return ENERGY_SCALE*dTim1;
  }

  // From here on, auxiliary functions
  vd Energy::kappaL() const {
    TRACE(0,"Energy::kappaL()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    
    // vd kappaL(v.gc->Ns,fillwith::zeros);
    // vd Tti=v.T.tdata();
    // vd kappait=v.gc->gas.kappa(Tti);

    // if(v.left==NULL){
    //   vd Ttip1=v.right->T.tdata();
    //   vd kappaitp1=v.gc->gas.kappa(Ttip1);
    //   kappaL=v.wL1*kappaitp1+v.wL0*kappait;
    // }
    // else{
    //   vd Ttim1=v.left->T.tdata();
    //   vd kappaitm1=v.gc->gas.kappa(Ttim1);
    //   kappaL=v.wLr*kappait+v.wLl*kappaitm1;
    // }
    // TRACE(100,"kappaL(0):"<<kappaL(0));
    // kappaL.zeros();		// WARNING
    vd kappaL(v.gc->Ns,fillwith::ones);
    kappaL*=v.gc->gas.kappa(v.gc->T0);
    return kappaL;
  }
  vd Energy::kappaR() const {		// Returns thermal conductivity time domain data
    TRACE(0,"Energy::kappaR()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    
    // vd kappaR(v.gc->Ns,fillwith::zeros);
    // vd Tti=v.T.tdata();
    // vd kappait=v.gc->gas.kappa(Tti);

    // if(v.right==NULL){
    //   vd Ttim1=v.left->T.tdata();
    //   vd kappaitm1=v.gc->gas.kappa(Ttim1);
    //   kappaR=v.wRNm2*kappaitm1+v.wRNm1*kappait;
    // }
    // else{
    //   vd Ttip1=v.right->T.tdata();
    //   vd kappaitp1=v.gc->gas.kappa(Ttip1);
    //   kappaR=v.wRl*kappait+v.wRr*kappaitp1;
    // }
    // kappaR.zeros();		// WARNING!!
    // TRACE(100,"kappaR(0):"<<kappaR(0));
    vd kappaR(v.gc->Ns,fillwith::ones);
    kappaR*=v.gc->gas.kappa(v.gc->T0);
    return kappaR;
  }
  d Energy::gamma() const {
    // d T0=v.T(0);
    d T0=v.gc->T0; 
    return v.gc->gas.gamma(T0);
  }
  Energy::~Energy(){
    TRACE(-5,"Energy destructor");
  }

} // namespace tube
