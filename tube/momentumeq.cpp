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
    cout << "--------- Showing momentum weight factors for i=" << vertex.i <<"\n" \
	 << "Wddt      : "<<Wddt      <<"\n"					\
	 << "Wuim1     : "<<Wuim1<<"\n"						\
	 << "Wui       : "<<Wui      <<"\n"					\
	 << "Wuip1     : "<<Wuip1<<"\n"						\
	 << "Wpim1     : "<<Wpim1      <<"\n"					\
	 << "Wpi       : "<<Wpi      <<"\n"					\
	 << "Wpip1     : "<<Wpip1      <<"\n"					\
            ;      
  }
  vd Momentum::Error(){		// Error in momentum equation
    TRACE(6,"Momentum::Error()");
    vd error(gc->Ns,fillwith::zeros);

    vd rhoti=vertex.rho.tdata();
    vd Uti=vertex.U.tdata();
    error+=Wddt*gc->DDTfd*gc->fDFT*(Uti%rhoti);
    error+=Wui*gc->fDFT*(rhoti%Uti%Uti);
    
    // Pressure terms    
    vd pi=vertex.p();
    error+=Wpi*pi;

    vd rhotim1(gc->Ns,fillwith::zeros);
    vd Utim1(gc->Ns,fillwith::zeros);
    vd pim1(gc->Ns,fillwith::zeros);
    vd Utip1(gc->Ns,fillwith::zeros);
    vd rhotip1(gc->Ns,fillwith::zeros);
    vd pip1(gc->Ns,fillwith::zeros);

    
    if(left!=NULL){
      rhotim1=left->rho.tdata();
      Utim1=left->U.tdata();
      error+=Wuim1*gc->fDFT*(rhotim1%Utim1%Utim1);
      // Pressure term
      pim1=left->p();
      error+=Wpim1*pim1;
    }
    if(right!=NULL){
      Utip1=right->U.tdata();
      rhotip1=right->rho.tdata();
      error+=Wuip1*gc->fDFT*(rhotip1%Utip1%Utip1);
      // Pressure term    
      pip1=right->p();
      error+=Wpip1*pip1;
    }

    // Artificial viscosity terms
    #ifdef MOM_VISCOSITY
    if(i>0 && i<Ncells-1){
      error+=-gc->rho0*D_r()*(right->U()-vertex.U());
      error+= gc->rho0*D_l()*(vertex.U()-left->U());
    }
    else if(i==0){
      error+=-gc->rho0*D_r()*(right->right->U()-right->U());
      error+= gc->rho0*D_l()*(right->U()-vertex.U());
    }
    else {
      error+=-gc->rho0*D_r()*(vertex.U()-left->U());
      error+= gc->rho0*D_l()*(left->U()-left->left->U());
    }
    #endif
    // Drag term
    // error+=vVf*tube.drag(i)/vSf;

    
    // (Boundary) source term
    error+=vertex.msource();
    error(0)*=MOM_SCALE0;
    return MOM_SCALE*error;
  }
  dmat Momentum::dUi(){
    TRACE(0,"Momentum::dUi()");
    dmat dUi=zero;
    // dUi+=vVf*tube.drag.dUi(i)/vSf;		       // Drag term
    dUi+=Wddt*gc->DDTfd*gc->fDFT*diagtmat(vertex.rho)*gc->iDFT; // Time-derivative term
    dUi+=2.0*Wui*gc->fDFT*(diagtmat(vertex.rho)*diagtmat(vertex.U))*gc->iDFT;
    // Artificial viscosity terms
    #ifdef MOM_VISCOSITY
    if(i>0 && i<Ncells-1){
      dUi+=gc->rho0*(D_l()+D_r());
    }
    else if(i==0)
      dUi+=-gc->rho0*D_l();
    else
      dUi+=-gc->rho0*D_r();
    #endif
    dUi.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dUi;
  }
  dmat Momentum::drhoi(){
    TRACE(0,"Momentum::drhoi()");
    dmat drhoi=zero;
    drhoi+=Wddt*gc->DDTfd*gc->fDFT*diagtmat(vertex.U)*gc->iDFT;
    drhoi+=Wui*gc->fDFT*diagtmat(vertex.U)*diagtmat(vertex.U)*gc->iDFT;
    drhoi.row(0)*=MOM_SCALE0;
    return MOM_SCALE*drhoi;
  }
  dmat Momentum::dpi(){
    TRACE(0,"Momentum::dpi()");
    dmat I(gc->Ns,gc->Ns,fillwith::eye);
    dmat dpi=zero;
    dpi+=Wpi*I;
    dpi.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dpi;
  }
  dmat Momentum::drhoim1(){
    TRACE(0,"Momentum::drhoim1()");
    dmat drhoim1=zero;
    if(left!=NULL)
      drhoim1+=Wuim1*gc->fDFT*diagtmat(left->U)*diagtmat(left->U)*gc->iDFT;
    drhoim1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*drhoim1;
  }
  dmat Momentum::dUim1(){
    TRACE(0,"Momentum::dUim1()");    // Todo: add this term!;
    dmat dUim1=zero;
    if(left!=NULL){
      dUim1+=2.0*Wuim1*gc->fDFT*diagtmat(left->rho)*diagtmat(left->U)*gc->iDFT;
    }
    #ifdef MOM_VISCOSITY
    // Artificial viscosity terms
    if(i>0 && i<Ncells-1){
      dUim1+=-gc->rho0*D_l();
    }
    else if(i==Ncells-1)
      dUim1+=gc->rho0*(D_r()+D_l());
    #endif
    
    dUim1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dUim1;
  }
  dmat Momentum::dpim1(){
    TRACE(0,"Momentum::dpim1()");
    dmat dpim1=zero;
    dmat I(gc->Ns,gc->Ns,fillwith::eye);
    if(left!=NULL)
      dpim1+=Wpim1*I;

    dpim1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dpim1;
  }
  dmat Momentum::drhoip1(){
    TRACE(0,"Momentum::dhoip1()");    // Todo: add this term!;
    dmat drhoip1=zero;
    if(right!=NULL)
      drhoip1+=Wuip1*gc->fDFT*diagtmat(right->U)*diagtmat(right->U)*gc->iDFT;

    drhoip1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*drhoip1;
  }
  dmat Momentum::dUip1(){
    TRACE(0,"Momentum::dUip1()"); // Todo: add this term!;
    dmat dUip1=zero;
    if(right!=NULL){
      dUip1+=2.0*Wuip1*gc->fDFT*diagtmat(right->rho)*diagtmat(right->U)*gc->iDFT;
    }
    #ifdef MOM_VISCOSITY
    // Artificial viscosity terms
    if(i>0 && i<Ncells-1){
      dUip1+=-gc->rho0*D_r();
    }
    else if(i==0)
      dUip1+=gc->rho0*(D_r()+D_l());
    #endif
    
    dUip1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dUip1;
  }
  dmat Momentum::dpip1(){
    TRACE(0,"Momentum::dpip1()");
    dmat dpip1=zero;
    dmat I(gc->Ns,gc->Ns,fillwith::eye);
    if(right!=NULL)
      dpip1+=Wpip1*I;

    dpip1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dpip1;
  }
  dmat Momentum::dUip2(){
    TRACE(0,"Momentum:dUip2()");
    // TRACE(50,"i:"<<i);
    dmat dUip2=zero;
    if(i==0 && left==NULL)
      // TRACE(20,"i is nul and left is nul");
      dUip2+=-gc->rho0*D_r();
    return dUip2;
  }
  dmat Momentum::dUim2(){
    dmat dUim2=zero;
    if((i==Ncells-1)&& right==NULL )
      // TRACE(20,"i is  Ncells and and right is nul");
      dUim2+=-gc->rho0*D_l();
    return dUim2;
  }  
  Momentum::~Momentum(){
    TRACE(-5,"Momentum destructor");
}
} // namespace tube
















