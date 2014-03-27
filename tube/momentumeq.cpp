#include "momentumeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{

  Momentum::Momentum(const Tube& tube,TubeVertex& gp):Equation(tube,gp){
    TRACE(0,"Momentum constructor...");
    // Standard boundary condition is an adiabatic no-slip wall
    if(i==0){			// Leftmost vertex
      TRACE(-1,"Leftmost vertex");
      Wuim1=0;
      Wui=wRl/SfR;
      Wuip1=wRr/SfR;
      Wpim1=0;
      Wpi=SfR*wRl-SfL*wL0+SfL-SfR;
      Wpip1=SfR*wRr-SfL*wL1;
      
    } else if(i==Ncells-1){	// Rightmost vertex
      TRACE(-1,"Rightmost vertex");
      Wuim1=-wLl/SfL;
      Wui=-wLr/SfL;
      Wuip1=0;

      Wpim1=-SfL*wLl+SfR*wRNm2;
      Wpi=-SfL*wLr+SfR*wRNm1+(SfL-SfR);
      Wpip1=0;

    } else{			// Normal interior vertex
      Wuim1=-wLl/SfL;
      Wui=wRl/SfR-wLr/SfL;
      Wuip1=wRr/SfR;

      Wpim1=-SfL*wLl;
      Wpi=SfL*(1.0-wLr)+SfR*(wRl-1.0);
      Wpip1=SfR*wRr;
    }
    TRACE(0,"Momentum constructor done");
  }
  vd Momentum::Error(){		// Error in momentum equation
    TRACE(0,"Momentum::Error()");
    vd error(Ns,fillwith::zeros);

    // WATCH IT! BELOW TERMS ARE ALL TIME DOMAIN DATA!! - EDIT: all but pressure is time domain data

    vd rhoti=vertex.rho.tdata();
    vd Uti=vertex.U.tdata();
    error+=vVf*DDTfd*(Uti%rhoti)/vSf;
    TRACE(-1,"Inbetween momentum error:"<< error);
    error+=Wui*fDFT*(rhoti%Uti%Uti);
    
    // Pressure terms    
    vd pip=vertex.p();	// pip because pi is already defined
    error+=Wpi*pip;

    if(i>0){
      vd rhotim1=left->rho.tdata();
      vd Utim1=left->U.tdata();
      error+=Wuim1*fDFT*(rhotim1%Utim1%Utim1);
      
      // Pressure terms   
      vd pim1=left->p();
      error+=Wpim1*pim1;

    }
    if(i<Ncells-1){
      vd Utip1=right->U.tdata();
      vd rhotip1=right->rho.tdata();
      error+=Wuip1*fDFT*(rhotip1%Utip1%Utip1);

      // Pressure terms    
      vd pip1=right->p();
      error+=Wpip1*pip1;
    }

    
    // Drag term
    error+=tube.drag(i);

    // (Boundary) source term
    error+=vertex.msource();
    return error;
  }
  dmat Momentum::dUi(){
    TRACE(0,"Momentum::dUi()");
    dmat dUi=zero;
    dUi+=tube.drag.dUi(i);		       // Drag term
    dUi+=DDTfd*fDFT*diagtmat(vertex.rho)*iDFT; // Time-derivative term
    dUi+=2.0*Wui*fDFT*(diagtmat(vertex.rho)*diagtmat(vertex.U))*iDFT;
    return dUi;
  }
  dmat Momentum::drhoi(){
    TRACE(0,"Momentum::drhoi()");
    dmat drhoi=zero;
    drhoi+=DDTfd*fDFT*diagtmat(vertex.U)*iDFT;
    drhoi+=Wui*fDFT*diagtmat(vertex.U)%diagtmat(vertex.U)*iDFT;
    return drhoi;
  }
  dmat Momentum::dpi(){
    TRACE(0,"Momentum::dpi()");
    dmat I(Ns,Ns,fillwith::eye);
    dmat dpi=zero;
    dpi+=Wpi*I;
    return dpi;
  }
  dmat Momentum::drhoim1(){
    TRACE(0,"Momentum::drhoim1()");
    dmat drhoim1=zero;
    if(i>0)
      drhoim1+=Wuim1*fDFT*diagtmat(left->U)*diagtmat(left->U)*iDFT;
    return drhoim1;
  }
  dmat Momentum::dUim1(){
    TRACE(0,"Momentum::dUim1()");    // Todo: add this term!;
    dmat dUim1=zero;
    if(i>0)
      dUim1+=2.0*Wuim1*fDFT*diagtmat(left->rho)*diagtmat(left->U)*iDFT;
    return dUim1;
  }
  dmat Momentum::dpim1(){
    TRACE(0,"Momentum::dpim1()");
    dmat dpim1=zero;
    dmat I(Ns,Ns,fillwith::eye);
    if(i>0)
      dpim1+=Wpim1*I;
    return dpim1;
  }
  dmat Momentum::drhoip1(){
    TRACE(0,"Momentum::dhoip1()");    // Todo: add this term!;
    dmat drhoip1=zero;
    if(i<Ncells-1)
      drhoip1+=Wuip1*fDFT*diagtmat(right->U)*diagtmat(right->U)*iDFT;
    return drhoip1;
  }
  dmat Momentum::dUip1(){
    TRACE(0,"Momentum::dUip1()"); // Todo: add this term!;
    dmat dUip1=zero;
    if(i<Ncells-1)
      dUip1+=2.0*Wuip1*fDFT*diagtmat(right->rho)*diagtmat(right->U)*iDFT;
    return dUip1;
  }
  dmat Momentum::dpip1(){
    TRACE(0,"Momentum::dpip1()");
    dmat dpip1=zero;
    dmat I(Ns,Ns,fillwith::eye);
    if(i<Ncells-1)
      dpip1+=Wpip1*I;
    return dpip1;
  }
  Momentum::~Momentum(){
    TRACE(-5,"Momentum destructor");
}
} // namespace tube





