#include "momentumeq.h"
#include "drag.h"
#include "momscale.h"
#include "tubevertex.h"
#define MOM_VISCOSITY

namespace tube{

  Momentum::Momentum(TubeVertex& gp):TubeEquation(gp){
    TRACE(0,"Momentum constructor...");
  }
  void Momentum::show(){
    cout << "--------- Showing momentum weight factors for i=" << v.i <<"\n" \
	 << "Wddt      : "<<Wddt      <<"\n"					\
	 << "Wuim1     : "<<Wuim1<<"\n"						\
	 << "Wui       : "<<Wui      <<"\n"					\
	 << "Wuip1     : "<<Wuip1<<"\n"						\
	 << "Wpim1     : "<<Wpim1      <<"\n"					\
	 << "Wpi       : "<<Wpi      <<"\n"					\
	 << "Wpip1     : "<<Wpip1      <<"\n"					\
            ;      
    v.gc->show();
  }
  vd Momentum::Error(){		// Error in momentum equation

    TRACE(6,"Momentum::Error()");
    vd error(v.gc->Ns,fillwith::zeros);

    vd rhoti=v.rho.tdata();
    vd Uti=v.U.tdata();
    error+=Wddt*v.gc->DDTfd*v.gc->fDFT*(Uti%rhoti);
    error+=Wui*v.gc->fDFT*(rhoti%Uti%Uti);
    
    // Pressure terms    
    vd pi=v.p();
    error+=Wpi*pi;

    vd rhotim1(v.gc->Ns,fillwith::zeros);
    vd Utim1(v.gc->Ns,fillwith::zeros);
    vd pim1(v.gc->Ns,fillwith::zeros);
    vd Utip1(v.gc->Ns,fillwith::zeros);
    vd rhotip1(v.gc->Ns,fillwith::zeros);
    vd pip1(v.gc->Ns,fillwith::zeros);

    
    if(v.left!=NULL){
      rhotim1=v.left->rho.tdata();
      Utim1=v.left->U.tdata();
      error+=Wuim1*v.gc->fDFT*(rhotim1%Utim1%Utim1);
      // Pressure term
      pim1=v.left->p();
      error+=Wpim1*pim1;
    }
    if(v.right!=NULL){
      Utip1=v.right->U.tdata();
      rhotip1=v.right->rho.tdata();
      error+=Wuip1*v.gc->fDFT*(rhotip1%Utip1%Utip1);
      // Pressure term    
      pip1=v.right->p();
      error+=Wpip1*pip1;
    }

    // Artificial viscosity terms
    #ifdef MOM_VISCOSITY
    if(v.i>0 && v.i<v.Ncells-1){
      error+=-v.gc->rho0*D_r()*(v.right->U()-v.U());
      error+= v.gc->rho0*D_l()*(v.U()-v.left->U());
    }
    else if(v.i==0){
      error+=-v.gc->rho0*D_r()*(v.right->right->U()-v.right->U());
      error+= v.gc->rho0*D_l()*(v.right->U()-v.U());
    }
    else {
      error+=-v.gc->rho0*D_r()*(v.U()-v.left->U());
      error+= v.gc->rho0*D_l()*(v.left->U()-v.left->left->U());
    }
    #endif
    // Drag term
    // error+=vVf*tube.drag(i)/vSf;

    
    // (Boundary) source term
    error+=v.msource();
    error(0)*=MOM_SCALE0;
    return MOM_SCALE*error;
  }
  dmat Momentum::dUi(){
    TRACE(0,"Momentum::dUi()");
    dmat dUi=zero;
    // dUi+=vVf*tube.drag.dUi(i)/vSf;		       // Drag term
    dUi+=Wddt*v.gc->DDTfd*v.gc->fDFT*diagtmat(v.rho)*v.gc->iDFT; // Time-derivative term
    dUi+=2.0*Wui*v.gc->fDFT*(diagtmat(v.rho)*diagtmat(v.U))*v.gc->iDFT;
    // Artificial viscosity terms
    #ifdef MOM_VISCOSITY
    if(v.i>0 && v.i<v.Ncells-1){
      dUi+=v.gc->rho0*(D_l()+D_r());
    }
    else if(v.i==0)
      dUi+=-v.gc->rho0*D_l();
    else
      dUi+=-v.gc->rho0*D_r();
    #endif
    dUi.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dUi;
  }
  dmat Momentum::drhoi(){
    TRACE(0,"Momentum::drhoi()");
    dmat drhoi=zero;
    drhoi+=Wddt*v.gc->DDTfd*v.gc->fDFT*diagtmat(v.U)*v.gc->iDFT;
    drhoi+=Wui*v.gc->fDFT*diagtmat(v.U)*diagtmat(v.U)*v.gc->iDFT;
    drhoi.row(0)*=MOM_SCALE0;
    return MOM_SCALE*drhoi;
  }
  dmat Momentum::dpi(){
    TRACE(0,"Momentum::dpi()");
    dmat I(v.gc->Ns,v.gc->Ns,fillwith::eye);
    dmat dpi=zero;
    dpi+=Wpi*I;
    dpi.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dpi;
  }
  dmat Momentum::drhoim1(){
    TRACE(0,"Momentum::drhoim1()");
    dmat drhoim1=zero;
    if(v.left!=NULL)
      drhoim1+=Wuim1*v.gc->fDFT*diagtmat(v.left->U)*diagtmat(v.left->U)*v.gc->iDFT;
    drhoim1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*drhoim1;
  }
  dmat Momentum::dUim1(){
    TRACE(0,"Momentum::dUim1()");    // Todo: add this term!;
    dmat dUim1=zero;
    if(v.left!=NULL){
      dUim1+=2.0*Wuim1*v.gc->fDFT*diagtmat(v.left->rho)*diagtmat(v.left->U)*v.gc->iDFT;
    }
    #ifdef MOM_VISCOSITY
    // Artificial viscosity terms
    if(v.i>0 && v.i<v.Ncells-1){
      dUim1+=-v.gc->rho0*D_l();
    }
    else if(v.i==v.Ncells-1)
      dUim1+=v.gc->rho0*(D_r()+D_l());
    #endif
    
    dUim1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dUim1;
  }
  dmat Momentum::dpim1(){
    TRACE(0,"Momentum::dpim1()");
    dmat dpim1=zero;
    dmat I(v.gc->Ns,v.gc->Ns,fillwith::eye);
    if(v.left!=NULL)
      dpim1+=Wpim1*I;

    dpim1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dpim1;
  }
  dmat Momentum::drhoip1(){
    TRACE(0,"Momentum::dhoip1()");    // Todo: add this term!;
    dmat drhoip1=zero;
    if(v.right!=NULL)
      drhoip1+=Wuip1*v.gc->fDFT*diagtmat(v.right->U)*diagtmat(v.right->U)*v.gc->iDFT;
    drhoip1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*drhoip1;
  }
  dmat Momentum::dUip1(){
    TRACE(0,"Momentum::dUip1()"); // Todo: add this term!;
    dmat dUip1=zero;
    if(v.right!=NULL){
      dUip1+=2.0*Wuip1*v.gc->fDFT*diagtmat(v.right->rho)*diagtmat(v.right->U)*v.gc->iDFT;
    }
    #ifdef MOM_VISCOSITY
    // Artificial viscosity terms
    if(v.i>0 && v.i<v.Ncells-1){
      dUip1+=-v.gc->rho0*D_r();
    }
    else if(v.i==0)
      dUip1+=v.gc->rho0*(D_r()+D_l());
    #endif
    
    dUip1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dUip1;
  }
  dmat Momentum::dpip1(){
    TRACE(0,"Momentum::dpip1()");
    dmat dpip1=zero;
    dmat I(v.gc->Ns,v.gc->Ns,fillwith::eye);
    if(v.right!=NULL)
      dpip1+=Wpip1*I;

    dpip1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dpip1;
  }
  dmat Momentum::dUip2(){
    TRACE(0,"Momentum:dUip2()");
    // TRACE(50,"i:"<<i);
    dmat dUip2=zero;
    if(v.i==0 && v.left==NULL)
      // TRACE(20,"i is nul and left is nul");
      dUip2+=-v.gc->rho0*D_r();
    return dUip2;
  }
  dmat Momentum::dUim2(){
    dmat dUim2=zero;
    if((v.i==v.Ncells-1)&& v.right==NULL )
      // TRACE(20,"i is  Ncells and and right is nul");
      dUim2+=-v.gc->rho0*D_l();
    return dUim2;
  }  
  Momentum::~Momentum(){
    TRACE(-5,"Momentum destructor");
}
} // namespace tube
















