#include "energyeq.h"
#include "tubevertex.h"
#include "conduction.h"

// #define ENERGY_SCALE (1.0/(v.gc->p0))
#define ENERGY_SCALE (1.0)
// #define EN_VISCOSITY


namespace tube{


  // void Energy::show() const{
  //   cout << "--------- Showing energy weight factors for i=" << v.i <<"\n" \
  // 	 << "v.eWddt      : "<<v.eWddt      <<"\n"					\
  // 	 << "v.eWgim1     : "<<v.eWgim1<<"\n"						\
  // 	 << "v.eWgi       : "<<v.eWgi      <<"\n"					\
  // 	 << "v.eWgip1     : "<<v.eWgip1<<"\n"						\
  // 	 << "v.eWjim1     : "<<v.eWjim1      <<"\n"					\
  // 	 << "v.eWji       : "<<v.eWji      <<"\n"					\
  // 	 << "v.eWjip1     : "<<v.eWjip1      <<"\n"					\
  // 	 << "v.eWc1     : "<<v.eWc1      <<"\n"					\
  // 	 << "v.eWc2     : "<<v.eWc2      <<"\n"					\
  // 	 << "v.eWc3     : "<<v.eWc3      <<"\n"					\
  // 	 << "v.eWc4     : "<<v.eWc4      <<"\n"	\

  //     ;      

  // }
  vd Energy::error(const TubeVertex& v) const {		// Error in momentum equation
    TRACE(6,"Energy::Error()");
    assert(v.gc!=NULL);
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    vd error(v.gc->Ns,fillwith::zeros);
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    vd Uti=v.U.tdata();
    vd pti=v.p.tdata()+v.getp0t();
    vd Tti=v.T.tdata();
    // error+=v.eWddt*DDTfd*v.p()/(gamma-1.0);
    // error+=v.eWgi*gamfac*fDFT*(pti%Uti);
    // error+=v.eWji*fDFT*(pti%Uti);
    error+=fDFT*((v.eWc2*kappaL(v)+v.eWc3*kappaR(v))%Tti);
    
    if(v.left!=NULL){
      // vd Utim1=v.left->U.tdata();
      // vd ptim1=v.left->p.tdata()+v.getp0t();
      // error+=v.eWgim1*gamfac*fDFT*(ptim1%Utim1);
      // error+=v.eWjim1*fDFT*(ptim1%Uti);
      // Conduction term
      vd Tim1=v.left->T.tdata();
      error+=v.eWc1*fDFT*(kappaL(v)%Tim1);

    }
    if(v.right!=NULL){
      // vd Utip1=v.right->U.tdata();
      // vd ptip1=v.right->p.tdata()+v.getp0t();
      // error+=v.eWgip1*gamfac*fDFT*(ptip1%Utip1);
      // error+=v.eWjip1*fDFT*(ptip1%Uti);
      vd Tip1=v.right->T.tdata();
      error+=v.eWc4*fDFT*(kappaR(v)%Tip1);
    }

    // Artificial viscosity terms      
    #ifdef EN_VISCOSITY
    const d& vSf=v.lg.vSf;
    if(v.i>0 && v.i<v.nCells-1){
      error+=-d_r(v)*(v.right->p() -v.p())*vSf;
      error+= d_l(v)*(v.p() -v.left->p())  *vSf;
    }
    else if(v.i==0){		// First v
      error+=-d_r(v)*(v.right->right->p()-v.right->p())*vSf;
      error+=d_l(v)*(v.right->p()-v.p())*vSf;
    }
    else {			// Last v
      error+=-d_r(v)*(v.p()-v.left->p())*vSf;
      error+=d_l(v)*(v.left->p()-v.left->left->p())*vSf;
    }
    #endif
    // (Boundary source term)
    error+=v.esource();
    return error;
  }
  dmat Energy::dpi(const TubeVertex& v) const {
    TRACE(0,"Energy::dpi()");
    d T0=v.T(0);
    d gamma=this->gamma(v);
    dmat dpi=v.zero;
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    // v.eWddt = vVf, as defined in tubevertex.cpp
    d gamfac=gamma/(gamma-1.0);
    dpi+=(v.eWddt/(gamma-1.0))*DDTfd;
    dpi+=v.eWgi*gamfac*fDFT*v.U.diagt()*iDFT;
    dpi+=v.eWji*fDFT*v.U.diagt()*iDFT;
    // Artificial viscosity terms    
    #ifdef EN_VISCOSITY
    const d& vVf=v.lg.vVf;
    const d& vSf=v.lg.vSf;
    if(v.i>0 && v.i<v.nCells-1){
      dpi+=(d_l(v)+d_r(v))*vSf;	// Middle vertex
    }
    else if(v.i==0)
      dpi+=-d_l(v)*vSf;	// First vertex
    else		
      dpi+=-d_r(v)*vSf;	// Last vertex
    #endif
    // return dpi;
    return v.zero;
  }
  dmat Energy::dpim1(const TubeVertex& v) const {
    TRACE(0,"Energy::dpim1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);    
    vd Uti=v.U.tdata();
    dmat dpim1=v.zero;
    if(v.left!=NULL){
      dpim1+=v.eWgim1*gamfac*fDFT*v.left->U.diagt()*iDFT;
    }
    dpim1+=v.eWjim1*fDFT*v.U.diagt()*iDFT;

    // Artificial viscosity terms
    #ifdef EN_VISCOSITY
    const d& vSf=v.lg.vSf;
    if(v.i>0 && v.i<v.nCells-1){
      dpim1+=-d_l(v)*vSf;
    }
    else if(v.i==v.nCells-1){		// Last vertex
      dpim1+=(d_l(v)+d_r(v))*vSf;
    }
    #endif
    // return dpim1;
    return v.zero;
  }
  dmat Energy::dpip1(const TubeVertex& v) const {
    TRACE(0,"Energy::dpip1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    vd Uti=v.U.tdata();
    dmat dpip1=v.zero;
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    if(v.right!=NULL){
      dpip1+=v.eWgip1*gamfac*fDFT*v.right->U.diagt()*iDFT;
    }
    dpip1+=v.eWjip1*fDFT*v.U.diagt()*iDFT;
    // Artificial viscosity terms
    #ifdef EN_VISCOSITY    
    const d& vSf=v.lg.vSf;
    if(v.i>0 && v.i<v.nCells-1){
      dpip1+=-d_r(v)*vSf;
    }
    else if(v.i==0){		// First vertex
      dpip1+=(d_l(v)+d_r(v))*vSf;
    }
    #endif
    // return dpip1;
    return v.zero;
  }
  dmat Energy::dUi(const TubeVertex& v) const {
    TRACE(0,"Energy::dUi()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    dmat dUi=v.zero;			    // Initialize with zeros
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    dmat diagpt=diagmat(v.getp0t()+v.p.tdata());
    dUi+=v.eWgi*gamfac*fDFT*diagpt*iDFT;
    dUi+=v.eWji*fDFT*diagpt*iDFT;
    if(v.left!=NULL)
      dUi+=v.eWjim1*fDFT*diagmat(v.left->p.tdata()+v.getp0t())*iDFT;
    if(v.right!=NULL)
      dUi+=v.eWjip1*fDFT*diagmat(v.right->p.tdata()+v.getp0t())*iDFT;
    // return dUi;
    return v.zero;
  }
  dmat Energy::dUim1(const TubeVertex& v) const {
    TRACE(0,"Energy::dUim1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    dmat dUim1=v.zero;
    if(v.left!=NULL){
      dUim1+=v.eWgim1*gamfac*fDFT*diagmat(v.getp0t()+v.left->p.tdata())*iDFT;
    }
    // return dUim1;
    return v.zero;
  }
  dmat Energy::dUip1(const TubeVertex& v) const {
    TRACE(0,"Energy::dUip1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    dmat dUip1=v.zero;
    if(v.right!=NULL){
      dUip1+=v.eWgip1*gamfac*fDFT*diagmat(v.getp0t()+v.right->p.tdata())*iDFT;
    }
    // return dUip1;
    return v.zero;
  }
  dmat Energy::dpip2(const TubeVertex& v) const {
    dmat dpip2=v.zero;
    #ifdef EN_VISCOSITY
    if(v.i==0 && v.left==NULL){
      const d& vSf=v.lg.vSf;
      dpip2+=-d_r(v)*vSf;
    }
    #endif 
    // return dpip2;
    return v.zero;
  }
  dmat Energy::dpim2(const TubeVertex& v) const {
    dmat dpim2=v.zero;
    #ifdef EN_VISCOSITY
    if(v.i==v.nCells-1 && v.right==NULL){
      const d& vSf=v.lg.vSf;
      dpim2+=-d_l(v)*v.lg.vSf;
    }
    #endif 
    // return dpim2; 
    return v.zero;   
  }
  // ############################## TEMPERATURE AND CONDUCTION TERMS
  dmat Energy::dTip1(const TubeVertex& v) const {
    TRACE(0,"Energy::dTip1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      

    dmat dTip1=v.zero;
    if(v.right!=NULL)    
      dTip1+=v.eWc4*fDFT*diagmat(kappaR(v))*iDFT;
    // dTip1.row(0)*=ENERGY_SCALE0;
    return dTip1;
  }
  dmat Energy::dTi(const TubeVertex& v) const {
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    TRACE(0,"Energy::dTi()");
    dmat dTi=v.zero;
    dTi+=fDFT*diagmat(v.eWc2*kappaL(v)+v.eWc3*kappaR(v))*iDFT;
    // dTi.row(0)*=ENERGY_SCALE0;
    return dTi;
  }
  dmat Energy::dTim1(const TubeVertex& v) const {
    TRACE(0,"Energy::dTim1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    dmat dTim1=v.zero;
    if(v.left!=NULL)    
      dTim1+=v.eWc1*fDFT*diagmat(kappaL(v))*iDFT;
    return dTim1;
  }

  // From here on, auxiliary functions
  vd Energy::kappaL(const TubeVertex& v)  const {
    TRACE(0,"Energy::kappaL(v)");
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
  vd Energy::kappaR(const TubeVertex& v)  const {		// Returns thermal conductivity time domain data
    TRACE(0,"Energy::kappaR(v)");
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
  d Energy::gamma(const TubeVertex& v) const {
    // d T0=v.T(0);
    d T0=v.gc->T0; 
    return v.gc->gas.gamma(T0);
  }


} // namespace tube
