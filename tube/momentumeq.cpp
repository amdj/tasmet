#include "momentumeq.h"
#include "tubevertex.h"
#include "tube.h"

#include "momscale.h"
#define MOM_VISCOSITY

namespace tube{

  Momentum::Momentum(const Tube& tube,TubeVertex& gp):TubeEquation(tube,gp){
    TRACE(0,"Momentum constructor...");
  }
  vd Momentum::Error(){		// Error in momentum equation
    TRACE(0,"Momentum::Error()");
    vd error(Ns,fillwith::zeros);

    vd rhoti=vertex.rho.tdata();
    vd Uti=vertex.U.tdata();
    error+=Wddt*DDTfd*fDFT*(Uti%rhoti);
    error+=Wui*fDFT*(rhoti%Uti%Uti);
    
    // Pressure terms    
    vd pi=vertex.p();
    error+=Wpi*pi;

    vd rhotim1(Ns,fillwith::zeros);
    vd Utim1(Ns,fillwith::zeros);
    vd pim1(Ns,fillwith::zeros);
    vd Utip1(Ns,fillwith::zeros);
    vd rhotip1(Ns,fillwith::zeros);
    vd pip1(Ns,fillwith::zeros);

    
    if(left!=NULL){
      rhotim1=left->rho.tdata();
      Utim1=left->U.tdata();
      error+=Wuim1*fDFT*(rhotim1%Utim1%Utim1);
      // Pressure term
      pim1=left->p();
      error+=Wpim1*pim1;
    }
    if(right!=NULL){
      Utip1=right->U.tdata();
      rhotip1=right->rho.tdata();
      error+=Wuip1*fDFT*(rhotip1%Utip1%Utip1);
      // Pressure term    
      pip1=right->p();
      error+=Wpip1*pip1;
    }

    // Artificial viscosity terms
    #ifdef MOM_VISCOSITY
    if(i>0 && i<Ncells-1){
      error+=-D_r()*(right->U()-vertex.U());
      error+= D_l()*(vertex.U()-left->U());
    }
    else if(i==0){
      error+=-D_r()*(right->right->U()-right->U());
      error+= D_l()*(right->U()-vertex.U());
    }
    else {
      error+=-D_r()*(vertex.U()-left->U());
      error+= D_l()*(left->U()-left->left->U());
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
    dUi+=Wddt*DDTfd*fDFT*diagtmat(vertex.rho)*iDFT; // Time-derivative term
    dUi+=2.0*Wui*fDFT*(diagtmat(vertex.rho)*diagtmat(vertex.U))*iDFT;

    // Artificial viscosity terms
    #ifdef MOM_VISCOSITY
    if(i>0 && i<Ncells-1){
      dUi+=(D_l()+D_r());
    }
    else if(i==0)
      dUi+=-D_l();
    else
      dUi+=-D_r();
    #endif
    
    dUi.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dUi;
  }
  dmat Momentum::drhoi(){
    TRACE(0,"Momentum::drhoi()");
    dmat drhoi=zero;
    drhoi+=Wddt*DDTfd*fDFT*diagtmat(vertex.U)*iDFT;
    drhoi+=Wui*fDFT*diagtmat(vertex.U)*diagtmat(vertex.U)*iDFT;
    drhoi.row(0)*=MOM_SCALE0;
    return MOM_SCALE*drhoi;
  }
  dmat Momentum::dpi(){
    TRACE(0,"Momentum::dpi()");
    dmat I(Ns,Ns,fillwith::eye);
    dmat dpi=zero;
    dpi+=Wpi*I;
    dpi.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dpi;
  }
  dmat Momentum::drhoim1(){
    TRACE(0,"Momentum::drhoim1()");
    dmat drhoim1=zero;
    if(left!=NULL)
      drhoim1+=Wuim1*fDFT*diagtmat(left->U)*diagtmat(left->U)*iDFT;
    drhoim1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*drhoim1;
  }
  dmat Momentum::dUim1(){
    TRACE(0,"Momentum::dUim1()");    // Todo: add this term!;
    dmat dUim1=zero;
    if(left!=NULL){
      dUim1+=2.0*Wuim1*fDFT*diagtmat(left->rho)*diagtmat(left->U)*iDFT;
    }
    #ifdef MOM_VISCOSITY
    // Artificial viscosity terms
    if(i>0 && i<Ncells-1){
      dUim1+=-D_l();
    }
    else if(i==Ncells-1)
      dUim1+=D_r()+D_l();
    #endif
    
    dUim1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dUim1;
  }
  dmat Momentum::dpim1(){
    TRACE(0,"Momentum::dpim1()");
    dmat dpim1=zero;
    dmat I(Ns,Ns,fillwith::eye);
    if(left!=NULL)
      dpim1+=Wpim1*I;

    dpim1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dpim1;
  }
  dmat Momentum::drhoip1(){
    TRACE(0,"Momentum::dhoip1()");    // Todo: add this term!;
    dmat drhoip1=zero;
    if(right!=NULL)
      drhoip1+=Wuip1*fDFT*diagtmat(right->U)*diagtmat(right->U)*iDFT;

    drhoip1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*drhoip1;
  }
  dmat Momentum::dUip1(){
    TRACE(0,"Momentum::dUip1()"); // Todo: add this term!;
    dmat dUip1=zero;
    if(right!=NULL){
      dUip1+=2.0*Wuip1*fDFT*diagtmat(right->rho)*diagtmat(right->U)*iDFT;
    }
    #ifdef MOM_VISCOSITY
    // Artificial viscosity terms
    if(i>0 && i<Ncells-1){
      dUip1+=-D_r();
    }
    else if(i==0)
      dUip1+=D_r()+D_l();
    #endif
    
    dUip1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dUip1;
  }
  dmat Momentum::dpip1(){
    TRACE(0,"Momentum::dpip1()");
    dmat dpip1=zero;
    dmat I(Ns,Ns,fillwith::eye);
    if(right!=NULL)
      dpip1+=Wpip1*I;

    dpip1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dpip1;
  }
  dmat Momentum::dUip2(){
    dmat dUip2=zero;
    if(i==0 && left==NULL)
      TRACE(20,"i is nul and left is nul");
      dUip2+=-D_r();
    return dUip2;
  }
  dmat Momentum::dUim2(){
    dmat dUim2=zero;
    if((i==Ncells-1)&& right==NULL )
      TRACE(20,"i is  Ncells and and right is nul");
      dUim2+=-D_l();
    return dUim2;
  }  
  Momentum::~Momentum(){
    TRACE(-5,"Momentum destructor");
}
} // namespace tube
















