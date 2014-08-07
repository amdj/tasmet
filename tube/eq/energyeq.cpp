#include "energyeq.h"
#include "tubevertex.h"


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
    TRACE(6,"Energy::Error(), i="<<v.i);
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
    vd rhoti=v.rho.tdata();

    
    error+=v.eWddt*DDTfd*v.p()/(gamma-1.0); // Static enthalpy term
    d vSfsq=pow(v.lg.vSf,2);
    error+=v.eWddt*DDTfd*fDFT*(0.5*rhoti%Uti%Uti/vSfsq); // This term should get a eWuddt factor later on
    error+=v.eWgi*fDFT*(gamfac*pti%Uti);
    // TRACE(100,"error:"<<error);
    error+=fDFT*(0.5*v.eWkini*rhoti%pow(Uti,3));
    error+=fDFT*((v.eWc2*kappaL(v)+v.eWc3*kappaR(v))%Tti);
    // TRACE(100,"error:"<<error);
    if(v.left!=NULL){
      vd Utim1=v.left->U.tdata();
      vd ptim1=v.left->p.tdata()+v.getp0t();
      vd rhotim1=v.left->rho.tdata();
      error+=fDFT*(v.eWgim1*gamfac*ptim1%Utim1);
      error+=fDFT*(v.eWkinim1*0.5*rhotim1%pow(Utim1,3));

      // Conduction term
      vd Tim1=v.left->T.tdata();
      error+=v.eWc1*fDFT*(kappaL(v)%Tim1);

    }
    // TRACE(100,"error:"<<error);
    if(v.right!=NULL){
      vd Utip1=v.right->U.tdata();
      vd ptip1=v.right->p.tdata()+v.getp0t();
      vd rhotip1=v.right->rho.tdata(); 
      error+=fDFT*(v.eWgip1*gamfac*ptip1%Utip1);
      error+=fDFT*(v.eWkinip1*0.5*rhotip1%pow(Utip1,3));
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
    // TRACE(100,"error:"<<error);
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
    return dpi;
  }
  dmat Energy::dpim1(const TubeVertex& v) const {
    TRACE(0,"Energy::dpim1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);    
    // vd Uti=v.U.tdata();
    dmat dpim1=v.zero;
    if(v.left!=NULL){
      dpim1+=v.eWgim1*gamfac*fDFT*v.left->U.diagt()*iDFT;
    }

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
    return dpim1;
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
    return dpip1;
  }
  dmat Energy::dUi(const TubeVertex& v) const {
    TRACE(0,"Energy::dUi()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    dmat dUi=v.zero;			    // Initialize with zeros
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);

    vd rhoti=v.rho.tdata();
    vd Uti=v.U.tdata();
    vd pti=v.p.tdata()+v.getp0t();

    d vSfsq=pow(v.lg.vSf,2);
    dUi+=DDTfd*fDFT*diagmat(v.eWddt*rhoti%Uti/vSfsq)*iDFT; // This term should get a eWuddt factor later on
    dUi+=fDFT*diagmat(v.eWgi*gamfac*pti)*iDFT;
    dUi+=fDFT*diagmat(v.eWkini*3*0.5*rhoti%pow(Uti,2))*iDFT;

    return dUi;
  }
  dmat Energy::drhoi(const TubeVertex& v) const {
    TRACE(0,"Energy::drhoi()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    dmat drhoi=v.zero;			    // Initialize with zeros
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);

    vd Uti=v.U.tdata();
    d vSfsq=pow(v.lg.vSf,2);
    drhoi+=DDTfd*fDFT*diagmat(0.5*v.eWddt*pow(Uti,2)/vSfsq)*iDFT; // This term should get a eWuddt factor later on
    drhoi+=fDFT*diagmat(v.eWkini*0.5*pow(Uti,3))*iDFT;

    return drhoi;
  }  
  dmat Energy::dUim1(const TubeVertex& v) const {
    TRACE(0,"Energy::dUim1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    dmat dUim1=v.zero;
    if(v.left!=NULL){
      vd Utim1=v.left->U.tdata();
      vd ptim1=v.left->p.tdata()+v.getp0t();
      vd rhotim1=v.left->rho.tdata(); 
      dUim1+=fDFT*diagmat(v.eWgim1*gamfac*ptim1)*iDFT;
      dUim1+=fDFT*diagmat(v.eWkinim1*3*0.5*rhotim1%pow(Utim1,2))*iDFT;

    }
    return dUim1;
  }
  dmat Energy::dUip1(const TubeVertex& v) const {
    TRACE(0,"Energy::dUip1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    dmat dUip1=v.zero;
    if(v.right!=NULL){
      vd Utip1=v.right->U.tdata();
      vd ptip1=v.right->p.tdata()+v.getp0t();
      vd rhotip1=v.right->rho.tdata(); 
      dUip1+=fDFT*diagmat(v.eWgip1*gamfac*ptip1)*iDFT;
      dUip1+=fDFT*diagmat(v.eWkinip1*3*0.5*rhotip1%pow(Utip1,2))*iDFT;

    }
    return dUip1;
  }
  dmat Energy::drhoim1(const TubeVertex& v) const {
    TRACE(0,"Energy::drhoim1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    dmat drhoim1=v.zero;
    if(v.left!=NULL){
      vd Utim1=v.left->U.tdata();
      drhoim1+=fDFT*diagmat(v.eWkinim1*0.5*pow(Utim1,3))*iDFT;
    }
    return drhoim1;
  }
  dmat Energy::drhoip1(const TubeVertex& v) const {
    TRACE(0,"Energy::drhoip1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    dmat drhoip1=v.zero;
    if(v.right!=NULL){
      vd rhotip1=v.right->rho.tdata();
      vd Utip1=v.right->U.tdata();
      drhoip1+=fDFT*diagmat(v.eWkinip1*0.5*pow(Utip1,3))*iDFT;
    }
    return drhoip1;
  }
  
  dmat Energy::dpip2(const TubeVertex& v) const {
    dmat dpip2=v.zero;
    #ifdef EN_VISCOSITY
    if(v.i==0 && v.left==NULL){
      const d& vSf=v.lg.vSf;
      dpip2+=-d_r(v)*vSf;
    }
    #endif 
    return dpip2;
  }
  dmat Energy::dpim2(const TubeVertex& v) const {
    dmat dpim2=v.zero;
    #ifdef EN_VISCOSITY
    if(v.i==v.nCells-1 && v.right==NULL){
      const d& vSf=v.lg.vSf;
      dpim2+=-d_l(v)*v.lg.vSf;
    }
    #endif 
    return dpim2; 
  }
  // ############################## TEMPERATURE AND CONDUCTION TERMS
  dmat Energy::dTip1(const TubeVertex& v) const {
    TRACE(0,"Energy::dTip1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      

    dmat dTip1=v.zero;
    if(v.right!=NULL)    
      {
	dTip1+=v.eWc4*fDFT*diagmat(kappaR(v))*iDFT;
      }
    return dTip1;
  }
  dmat Energy::dTi(const TubeVertex& v) const {
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    TRACE(0,"Energy::dTi()");
    dmat dTi=v.zero;
    dTi+=fDFT*diagmat(v.eWc2*kappaL(v)+v.eWc3*kappaR(v))*iDFT;
    return dTi;
  }
  dmat Energy::dTim1(const TubeVertex& v) const {
    TRACE(0,"Energy::dTim1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    dmat dTim1=v.zero;
    if(v.left!=NULL){
      dTim1+=v.eWc1*fDFT*diagmat(kappaL(v))*iDFT;
    }
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
