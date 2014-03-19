#include "momentumeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{

  Momentum::Momentum(const Tube& tube,const TubeVertex& gp):Equation(tube,gp){
    TRACE(0,"Momentum constructor done");
  }
  dmat Momentum::operator()(){
    TRACE(0,"Momentum::operator()");
    return Equation::operator()();
  }
  vd Momentum::Error(){		// Error in momentum equation
    TRACE(0,"Momentum::Error()");
    vd error(Ns,fillwith::zeros);

    // WATCH IT! ABOVE TERMS ARE ALL TIME DOMAIN DATA!!
    vd rhoi=vertex.rho.tdata();
    vd Ui=vertex.U.tdata();
    vd rhoip1=tube.vvertex[i+1].rho.tdata();
    vd rhoim1=tube.vvertex[i-1].rho.tdata();
    vd Uip1=tube.vvertex[i+1].U.tdata();
    vd Uim1=tube.vvertex[i-1].U.tdata();
    vd pim1=tube.vvertex[i-1].p.tdata();
    vd pip1=tube.vvertex[i+1].p.tdata();
    vd pit=vertex.p.tdata();	// pit because pi is already defined

    error+=vVf*DDTfd*(vertex.U*vertex.rho).tdata()/vSf;
    error+=(wRr/vSf)*fDFT*(rhoip1%Uip1%Uip1);
    error+=(wRl/SfR-wLr/SfL)*fDFT*(rhoi%Ui%Ui);
    error+=-1.0*(wLl/SfL)*fDFT*(rhoim1%Uim1%Uim1);

    // Pressure terms
    error+=SfR*wRr*pip1;
    error+=(SfR*(wRl+1.0)-SfL*(wLr+1.0))*pit;
    error+=-1.0*SfL*wLl*pim1;
    
    // Drag term
    error+=tube.drag(i);
    return error;
  }
  dmat Momentum::dUi(){
    TRACE(0,"Momentum::dUi()");
    d fac1=(wRl/SfR-wLr/SfL);
    dmat ddragdU=tube.drag.dUi(i);
    return ddragdU+DDTfd*fDFT*diagtmat(vertex.rho)*iDFT+fac1*fDFT*(diagtmat(vertex.rho)*diagtmat(vertex.U));
  }
  dmat Momentum::drhoi(){
    TRACE(0,"Momentum::drhoi()");
    d fac1=(wRl/SfR-wLr/SfL);
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





