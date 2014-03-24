#include "momentumeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{

  Momentum::Momentum(const Tube& tube,const TubeVertex& gp):Equation(tube,gp){
    TRACE(0,"Momentum constructor done");
    // Standard boundary condition is an adiabatic no-slip wall
    if(i==0){			// Leftmost vertex
      
    } else if(i==Ncells-1){	// Rightmost vertex

    } else{			// Normal interior vertex

    }
  }
  dmat Momentum::operator()(){
    TRACE(0,"Momentum::operator()");
    return Equation::operator()();
  }
  vd Momentum::Error(){		// Error in momentum equation
    TRACE(0,"Momentum::Error()");
    vd error(Ns,fillwith::zeros);

    // WATCH IT! BELOW TERMS ARE ALL TIME DOMAIN DATA!!
    vd rhoti=vertex.rho.tdata();
    vd Uti=vertex.U.tdata();
    vd rhotip1=tube.vvertex[i+1].rho.tdata();
    vd rhotim1=tube.vvertex[i-1].rho.tdata();
    vd Uitp1=tube.vvertex[i+1].U.tdata();
    vd Uitm1=tube.vvertex[i-1].U.tdata();

    error+=vVf*DDTfd*(vertex.U*vertex.rho).tdata()/vSf;
    error+=Wuip1*fDFT*(rhoip1%Uip1%Uip1);
    error+=Wui*fDFT*(rhoi%Ui%Ui);
    error+=Wuim1*fDFT*(rhoim1%Uim1%Uim1);

    // Pressure terms    
    // Frequency domain data
    vd pim1=tube.vvertex[i-1].p();
    vd pip1=tube.vvertex[i+1].p.tdata();
    vd pip=vertex.p();	// pip because pi is already defined

    error+=Wpip1*pip1;
    error+=Wpi*pip;
    error+=Wpim1*pim1;
    
    // Drag term
    error+=tube.drag(i);
    return error;
  }
  dmat Momentum::dUi(){
    TRACE(0,"Momentum::dUi()");
    dmat dUi=zero;
    dmat dUi+=tube.drag.dUi(i);		       // Drag term
    dUi+=DDTfd*fDFT*diagtmat(vertex.rho)*iDFT; // Time-derivative term
    dUi+=2.0*Wui*fDFT*(diagtmat(vertex.rho)*diagtmat(vertex.U))*iDFT;
    return dUi;
  }
  dmat Momentum::drhoi(){
    TRACE(0,"Momentum::drhoi()");
    dmat drhoi=zero;
    drhoi+=DDTfd*fDFT*Uid*iDFT;
    drhoi+=Wui*fDFT*diagtmat(vertex.U)%diagtmat(vertex.U)*iDFT;
    return drhoi;
  }
  dmat Momentum::dpi(){
    TRACE(0,"Momentum::dpi()");
    dmat I(Ns,Ns,fillwith::eye);
    dmat dpi=zero;
    dpi+=wpi*I;
    return dpi;
  }
  dmat Momentum::drhoim1(){
    TRACE(0,"Momentum::drhoim1()");
    dmat drhoim1=zero;
    if(i>0)
      drhoim1+=Wuim1*fDFT*diagtmat(tube.vvertex[i-1].U)*diagtmat(tube.vvertex[i-1].U)*iDFT;
    return drhoim1;
  }
  dmat Momentum::dUim1(){
    TRACE(0,"Momentum::dUim1()");    // Todo: add this term!;
    dmat dUim1=zero;
    if(i>0)
      dUim1+=-2.0*Wuim1*fDFT*diagtmat(tube.vvertex[i-1].rho)*diagtmat(tube.vvertex[i-1].U)*iDFT;
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
      drhoip1+=Wuip1*fDFT*diagtmat(tube.vvertex[i+1].U)*diagtmat(tube.vvertex[i+1].U)*iDFT;
    return drhoip1;
  }
  dmat Momentum::dUip1(){
    TRACE(0,"Momentum::dUip1()"); // Todo: add this term!;
    dmat dUip1=zero;
    if(i<Ncells-1)
      2.0*fDFT*diagtmat(tube.vvertex[i+1].rho)*diagtmat(tube.vvertex[i+1].U)*iDFT;
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





