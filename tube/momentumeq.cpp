#include "momentumeq.h"
#include "tubevertex.h"
#include "tube.h"

#include "momscale.h"
#define MOM_VISCOSITY

namespace tube{

  Momentum::Momentum(const Tube& tube,TubeVertex& gp):TubeEquation(tube,gp){
    TRACE(0,"Momentum constructor...");
    // Standard boundary condition is an adiabatic no-slip wall


    if(i==0){			// Leftmost vertex
      Wuim1=0;
      Wui=wRl/SfR;
      Wuip1=wRr/SfR;
      Wpim1=0;
      Wpi=SfR*wRl-SfL*wL0;
      Wpip1=SfR*wRr-SfL*wL1;
      
    } else if(i==Ncells-1){	// Rightmost vertex
      Wuim1=-wLl/SfL;
      Wui=-wLr/SfL;
      Wuip1=0;

      Wpim1=-SfL*wLl+SfR*wRNm2;
      Wpi=-SfL*wLr+SfR*wRNm1;
      Wpip1=0;

      
    } else{			// Normal interior vertex
      Wuim1=-wLl/SfL;
      Wui=wRl/SfR-wLr/SfL;
      Wuip1=wRr/SfR;

      Wpim1=-SfL*wLl;
      Wpi  = SfR*wRl-SfL*wLr;
      Wpip1=SfR*wRr;
      
    }
    // Contribution from changing cross-sectional area
    Wpi+=SfL-SfR;
    TRACE(0,"Momentum constructor done");
  }
  vd Momentum::Error(){		// Error in momentum equation
    TRACE(0,"Momentum::Error()");
    vd error(Ns,fillwith::zeros);

    vd rhoti=vertex.rho.tdata();
    vd Uti=vertex.U.tdata();
    error+=vVf*DDTfd*fDFT*(Uti%rhoti)/vSf;
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

    
    if(i>0){
      rhotim1=left->rho.tdata();
      Utim1=left->U.tdata();
      error+=Wuim1*fDFT*(rhotim1%Utim1%Utim1);
      // Pressure term
      pim1=left->p();
      error+=Wpim1*pim1;
    }
    if(i<Ncells-1){
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
    dUi+=vVf*DDTfd*fDFT*diagtmat(vertex.rho)*iDFT/vSf; // Time-derivative term
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
    drhoi+=vVf*DDTfd*fDFT*diagtmat(vertex.U)*iDFT/vSf;
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
    if(i>0)
      drhoim1+=Wuim1*fDFT*diagtmat(left->U)*diagtmat(left->U)*iDFT;
    drhoim1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*drhoim1;
  }
  dmat Momentum::dUim1(){
    TRACE(0,"Momentum::dUim1()");    // Todo: add this term!;
    dmat dUim1=zero;
    if(i>0){
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
    if(i>0)
      dpim1+=Wpim1*I;

    dpim1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dpim1;
  }
  dmat Momentum::drhoip1(){
    TRACE(0,"Momentum::dhoip1()");    // Todo: add this term!;
    dmat drhoip1=zero;
    if(i<Ncells-1)
      drhoip1+=Wuip1*fDFT*diagtmat(right->U)*diagtmat(right->U)*iDFT;

    drhoip1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*drhoip1;
  }
  dmat Momentum::dUip1(){
    TRACE(0,"Momentum::dUip1()"); // Todo: add this term!;
    dmat dUip1=zero;
    if(i<Ncells-1){
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
    if(i<Ncells-1)
      dpip1+=Wpip1*I;

    dpip1.row(0)*=MOM_SCALE0;
    return MOM_SCALE*dpip1;
  }
  dmat Momentum::dUip2(){
    dmat dUip2=zero;
    if(i==0)
      dUip2+=-D_r();
    return dUip2;
  }
  dmat Momentum::dUim2(){
    dmat dUim2=zero;
    if(i==Ncells-1)
      dUim2+=-D_l();
    return dUim2;
  }  
  Momentum::~Momentum(){
    TRACE(-5,"Momentum destructor");
}
} // namespace tube
















