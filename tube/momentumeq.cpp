#include "momentumeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{

  Momentum::Momentum(const Tube& tube,const TubeVertex& gp):Equation(tube,gp){
    TRACE(0,"Momentum constructor done");
    // Standard boundary condition is an adiabatic no-slip wall
    if(i==0){
      
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
------------------------IK BEN HIER
    dmat Uid=diagtmat(vertex.U);
    return DDTfd*fDFT*Uid*iDFT+fac1*fDFT*Uid*Uid*iDFT;
  }
  dmat Momentum::dpi(){
    TRACE(0,"Momentum::dpi()");
    dmat I(Ns,Ns,fillwith::eye);
    return (SfR*(wRl+1.0)-SfL*(wLr+1))*I;
  }
  dmat Momentum::drhoim1(){
    TRACE(0,"Momentum::drhoim1()");
    return -1.0*(wLl/SfL)*fDFT*diagtmat(tube.vvertex[i-1].U)*diagtmat(tube.vvertex[i-1].U)*iDFT;
  }
  dmat Momentum::dUim1(){
    TRACE(0,"Momentum::dUim1()")    // Todo: add this term!;
    return -2.0*(wLl/SfL)*fDFT*diagtmat(tube.vvertex[i-1].rho)*diagtmat(tube.vvertex[i-1].U)*iDFT;
  }
  dmat Momentum::dpim1(){
    TRACE(0,"Momentum::dpim1()");;
    dmat I(Ns,Ns,fillwith::eye);
    return -1.0*(wLl*SfL)*I;
  }
  dmat Momentum::drhoip1(){
    TRACE(0,"Momentum::drhoip1()")    // Todo: add this term!;
    return 1.0*(wRr/SfR)*fDFT*diagtmat(tube.vvertex[i+1].U)*diagtmat(tube.vvertex[i+1].U)*iDFT;
  }
  dmat Momentum::dUip1(){
    TRACE(0,"Momentum::dUip1()")    // Todo: add this term!;
    return 2.0*(wRr/SfR)*fDFT*diagtmat(tube.vvertex[i+1].rho)*diagtmat(tube.vvertex[i+1].U)*iDFT;
  }
  dmat Momentum::dpip1(){
    TRACE(0,"Momentum::dpip1()");
    dmat I(Ns,Ns,fillwith::eye);
    return 1.0*(wRr*SfR)*I;
  }
  Momentum::~Momentum(){
    TRACE(-5,"Momentum destructor");
}
} // namespace tube





