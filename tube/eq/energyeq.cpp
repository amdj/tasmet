#include "energyeq.h"
#include "tubevertex.h"
#include "tube.h"

#include "artvisco.h"


namespace tube{

  void Energy::show() const{
    cout << "Full energy equation\n";
    #ifdef EN_VISCOSITY
    cout << "Artificial viscosity turned ON for energy equation\n";
    #else
    cout << "Artificial viscosity turned OFF for energy equation\n";
    #endif
  }
  void Energy::init(const Tube& t){
    TRACE(8,"Energy::init(tube)");
    heat=&t.getHeatSource();
  }
  JacRow Energy::jac(const TubeVertex& v) const{
    TRACE(6,"Energy::jac()");
    JacRow jac(dofnr,12);
    jac+=drhoi(v);
    jac+=dTi(v);
    jac+=dpi(v);
    jac+=dUi(v);

    if(v.left){
      jac+=drhoim1(v);
      jac+=dTim1(v);
      jac+=dpim1(v);
      jac+=dUim1(v);
    }
    if(v.right){
      jac+=drhoip1(v);
      jac+=dTip1(v);
      jac+=dpip1(v);
      jac+=dUip1(v);
    }
    return jac;    
  }
  vd Energy::error(const TubeVertex& v) const {		// Error in momentum equation
    TRACE(6,"Energy::Error(), i="<<v.i);
    assert(v.gc!=NULL);
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    vd error(v.gc->Ns,fillwith::zeros);
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);

    const vd& Uti=v.U.tdata();
    vd pti=0.5*(v.pL().tdata()+v.pR().tdata())+v.getp0t();
    vd pRt=v.getp0t()+v.pR();
    vd pLt=v.getp0t()+v.pL();
    const vd& Tti=v.T.tdata();
    const vd& rhoti=v.rho.tdata();

    error+=v.eWddt*DDTfd*0.5*(v.pL()()+v.pR()())/(gamma-1.0); // Static
    // enthalpy term
    
    error+=v.eWddtkin*DDTfd*fDFT*(rhoti%Uti%Uti); // This term should get a eWddtkin factor later on
    error+=v.eWgi*fDFT*(gamfac*pti%Uti);
    error+=fDFT*(v.eWkini*rhoti%pow(Uti,3));
    error+=fDFT*((v.eWc2*kappaL(v)+v.eWc3*kappaR(v))%Tti);

    if(v.left){
      const vd& Utim1=v.left->U.tdata();
      vd ptim1=v.left->p.tdata()+v.getp0t();
      const vd& rhotim1=v.left->rho.tdata();
      error+=fDFT*(v.eWgim1*gamfac*ptim1%Utim1);
      error+=fDFT*(v.eWkinim1*rhotim1%pow(Utim1,3));

      // Conduction term
      const vd& Tim1=v.left->T.tdata();
      error+=v.eWc1*fDFT*(kappaL(v)%Tim1);

    }
    // TRACE(100,"error:"<<error);
    if(v.right){
      const vd& Utip1=v.right->U.tdata();
      vd ptip1=v.right->p.tdata()+v.getp0t();
      const vd& rhotip1=v.right->rho.tdata(); 
      error+=fDFT*(v.eWgip1*gamfac*ptip1%Utip1);
      error+=fDFT*(v.eWkinip1*rhotip1%pow(Utip1,3));
      const vd& Tip1=v.right->T.tdata();
      error+=v.eWc4*fDFT*(kappaR(v)%Tip1);
    }

    // Artificial viscosity terms      
    #ifdef EN_VISCOSITY
    TRACE(10,"Energy viscosity defined");
    if(v.left!=NULL && v.right!=NULL){
      error+=d_l(v)*(v.cWart2*v.p()+v.cWart1*v.left->p() );
      error+=d_r(v)*(v.cWart3*v.p()+v.cWart4*v.right->p() );
    }
    #endif
    assert(heat!=NULL);
    error+=v.eWddt*heat->heat(v);
    // (Boundary source term)
    error+=v.esource();
    // TRACE(100,"error:"<<error);
    return error;
  }
  vd Energy::domg(const TubeVertex& v) const {
    TRACE(0,"Energy::domg()");
    assert(v.gc!=NULL);
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);

    vd domg=v.eWddt*DDTfd*v.p()/(gamma-1.0)/v.gc->getomg(); // Static enthalpy term
    domg+=v.eWddtkin*DDTfd*fDFT*(v.rho.tdata()%v.U.tdata()%v.U.tdata())/v.gc->getomg();
    return domg;
  }

  JacCol Energy::dpi(const TubeVertex& v) const {
    TRACE(0,"Energy::dpi()");
    d T0=v.T(0);
    d gamma=this->gamma(v);
    JacCol dpi(v.p);
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    // v.eWddt = vVf, as defined in tubevertex.cpp
    d gamfac=gamma/(gamma-1.0);
    dpi+=(v.eWddt/(gamma-1.0))*DDTfd;
    dpi+=v.eWgi*gamfac*fDFT*v.U.diagt()*iDFT;

    // Artificial viscosity terms    
    #ifdef EN_VISCOSITY
    if(v.left!=NULL && v.right!=NULL){
      dpi+=d_l(v)*v.cWart2+d_r(v)*v.cWart3;	// Middle vertex
    }
    #endif
    return dpi;
  }
  JacCol Energy::dpim1(const TubeVertex& v) const {
    TRACE(0,"Energy::dpim1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);    
    // vd Uti=v.U.tdata();
    JacCol dpim1(v.left->p);
    if(v.left!=NULL){
      dpim1+=v.eWgim1*gamfac*fDFT*v.left->U.diagt()*iDFT;
    }

    // Artificial viscosity terms
    #ifdef EN_VISCOSITY
    if(v.left!=NULL && v.right!=NULL){
      dpim1+=d_l(v)*v.cWart1;
    }
    #endif
    return dpim1;
  }
  JacCol Energy::dpip1(const TubeVertex& v) const {
    TRACE(0,"Energy::dpip1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    const vd& Uti=v.U.tdata();
    JacCol dpip1(v.right->p);
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    dpip1+=v.eWgip1*gamfac*fDFT*v.right->U.diagt()*iDFT;

    // Artificial viscosity terms
    #ifdef EN_VISCOSITY    
    if(v.left!=NULL && v.right!=NULL){
      dpip1+=d_r(v)*v.cWart4;
    }
    #endif
    return dpip1;
  }
  JacCol Energy::dUi(const TubeVertex& v) const {
    TRACE(0,"Energy::dUi()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    JacCol dUi(v.U);			    // Initialize with zeros
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);

    const vd& rhoti=v.rho.tdata();
    const vd& Uti=v.U.tdata();
    vd pti=v.p.tdata()+v.getp0t();
    dUi+=DDTfd*fDFT*diagmat(2*v.eWddtkin*rhoti%Uti)*iDFT; // This term should get a eWuddt factor later on
    dUi+=fDFT*diagmat(v.eWgi*gamfac*pti)*iDFT;
    dUi+=fDFT*diagmat(v.eWkini*3*rhoti%pow(Uti,2))*iDFT;
    assert(heat!=NULL);
    dUi+=v.eWddt*heat->dUi(v);
    return dUi;
  }
  JacCol Energy::drhoi(const TubeVertex& v) const {
    TRACE(0,"Energy::drhoi()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    JacCol drhoi(v.rho);			    // Initialize with zeros
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);

    const vd& Uti=v.U.tdata();
    drhoi+=DDTfd*fDFT*diagmat(v.eWddtkin*pow(Uti,2))*iDFT; // This term should get a eWuddt factor later on
    drhoi+=fDFT*diagmat(v.eWkini*pow(Uti,3))*iDFT;

    return drhoi;
  }  
  JacCol Energy::dUim1(const TubeVertex& v) const {
    TRACE(0,"Energy::dUim1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    JacCol dUim1(v.left->U);
    const vd& Utim1=v.left->U.tdata();
    vd ptim1=v.left->p.tdata()+v.getp0t();
    const vd& rhotim1=v.left->rho.tdata(); 
    dUim1+=fDFT*diagmat(v.eWgim1*gamfac*ptim1)*iDFT;
    dUim1+=fDFT*diagmat(v.eWkinim1*3*rhotim1%pow(Utim1,2))*iDFT;

    return dUim1;
  }
  JacCol Energy::dUip1(const TubeVertex& v) const {
    TRACE(0,"Energy::dUip1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    JacCol dUip1(v.right->U);
    const vd& Utip1=v.right->U.tdata();
    vd ptip1=v.right->p.tdata()+v.getp0t();
    const vd& rhotip1=v.right->rho.tdata(); 
    dUip1+=fDFT*diagmat(v.eWgip1*gamfac*ptip1)*iDFT;
    dUip1+=fDFT*diagmat(v.eWkinip1*3*rhotip1%pow(Utip1,2))*iDFT;

    return dUip1;
  }
  JacCol Energy::drhoim1(const TubeVertex& v) const {
    TRACE(0,"Energy::drhoim1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    JacCol drhoim1(v.left->rho);
    const vd& Utim1=v.left->U.tdata();
    drhoim1+=fDFT*diagmat(v.eWkinim1*pow(Utim1,3))*iDFT;
    return drhoim1;
  }
  JacCol Energy::drhoip1(const TubeVertex& v) const {
    TRACE(0,"Energy::drhoip1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    d gamma=this->gamma(v);
    d gamfac=gamma/(gamma-1.0);
    JacCol drhoip1(v.right->rho);

    const vd& rhotip1=v.right->rho.tdata();
    const vd& Utip1=v.right->U.tdata();
    drhoip1+=fDFT*diagmat(v.eWkinip1*pow(Utip1,3))*iDFT;

    return drhoip1;
  }
  
  // JacCol Energy::dpip2(const TubeVertex& v) const {
  //   dmat dpip2=v.zero;
  //   #ifdef EN_VISCOSITY
  //   // if(v.i==0 && v.left==NULL){
  //   //   const d& vSf=v.lg.vSf;
  //   //   dpip2+=-d_r(v)*vSf;
  //   // }
  //   #endif 
  //   return dpip2;
  // }
  // JacCol Energy::dpim2(const TubeVertex& v) const {
  //   dmat dpim2=v.zero;
  //   #ifdef EN_VISCOSITY
  //   // if(v.i==v.nCells-1 && v.right==NULL){
  //   //   const d& vSf=v.lg.vSf;
  //   //   dpim2+=-d_l(v)*v.lg.vSf;
  //   // }
  //   #endif 
  //   return dpim2; 
  // }
  // ############################## TEMPERATURE AND CONDUCTION TERMS
  JacCol Energy::dTip1(const TubeVertex& v) const {
    TRACE(0,"Energy::dTip1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      

    JacCol dTip1(v.right->T);
    // TRACE(100,"dTip1, right is"<<v.right);
    dTip1+=v.eWc4*fDFT*diagmat(kappaR(v))*iDFT;

    return dTip1;
  }
  JacCol Energy::dTi(const TubeVertex& v) const {
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    TRACE(0,"Energy::dTi()");
    JacCol dTi(v.T);
    dTi+=fDFT*diagmat(v.eWc2*kappaL(v)+v.eWc3*kappaR(v))*iDFT;
    assert(heat!=NULL);
    dTi+=v.eWddt*heat->dTi(v);
    return dTi;
  }
  JacCol Energy::dTim1(const TubeVertex& v) const {
    TRACE(0,"Energy::dTim1()");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    JacCol dTim1(v.left->T);
    dTim1+=v.eWc1*fDFT*diagmat(kappaL(v))*iDFT;
    return dTim1;
  }

  // From here on, auxiliary functions
  vd Energy::kappaL(const TubeVertex& v)  const {
    TRACE(0,"Energy::kappaL(v)");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    
    vd kappaL(v.gc->Ns,fillwith::zeros);
    vd Tti=v.T.tdata();
    vd kappait=v.gc->gas.kappa(Tti);

    if(v.left==NULL){
      vd Ttip1=v.right->T.tdata();
      vd kappaitp1=v.gc->gas.kappa(Ttip1);
      kappaL=v.w.wL1*kappaitp1+v.w.wL0*kappait;
    }
    else{
      const vd& Ttim1=v.left->T.tdata();
      vd kappaitm1=v.gc->gas.kappa(Ttim1);
      kappaL=v.w.wLr*kappait+v.w.wLl*kappaitm1;
    }
    // TRACE(100,"kappaL(0):"<<kappaL(0));
    // kappaL.zeros();		// WARNING
    // vd kappaL(v.gc->Ns,fillwith::ones);
    // kappaL*=v.gc->gas.kappa(v.gc->T0);
    return kappaL;
  }
  vd Energy::kappaR(const TubeVertex& v)  const {		// Returns thermal conductivity time domain data
    TRACE(0,"Energy::kappaR(v)");
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;      
    
    vd kappaR(v.gc->Ns,fillwith::zeros);
    const vd& Tti=v.T.tdata();
    vd kappait=v.gc->gas.kappa(Tti);

    if(v.right==NULL){
      const vd& Ttim1=v.left->T.tdata();
      vd kappaitm1=v.gc->gas.kappa(Ttim1);
      kappaR=v.w.wRNm2*kappaitm1+v.w.wRNm1*kappait;
    }
    else{
      const vd& Ttip1=v.right->T.tdata();
      vd kappaitp1=v.gc->gas.kappa(Ttip1);
      kappaR=v.w.wRl*kappait+v.w.wRr*kappaitp1;
    }
    // kappaR.zeros();		// WARNING!!
    // TRACE(100,"kappaR(0):"<<kappaR(0));
    // vd kappaR(v.gc->Ns,fillwith::ones);
    // kappaR*=v.gc->gas.kappa(v.gc->T0);
    return kappaR;
  }
  d Energy::gamma(const TubeVertex& v) const {
    // d T0=v.T(0);
    d T0=v.gc->T0; 
    return v.gc->gas.gamma(T0);
  }


} // namespace tube
