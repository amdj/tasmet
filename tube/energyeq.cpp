#include "energyeq.h"
#include "vertex.h"
#include "tube.h"

namespace tube{

  Energy::Energy(const Tube& tube,TubeVertex& gp):Equation(tube,gp){
    // Standard boundary condition is adiabatic-no-slip-wall
    if(i==0){			// Leftmost vertex
      TRACE(-1,"Leftmost vertex");
      d xip1=tube.geom.vx(i+1);
      d xi=tube.geom.vx(i);
      d dxp=xip1-xi;

      Whim1=0;
      Whi=wRl-wL0;
      Whip1=wRr-wL1;

      Wjim1=0;
      Wji=wL0-wRl;
      Wjip1=wL1-wRr;
      
      Wc1=0;
      Wc2=0;
      Wc3=SfR/dxp;
      Wc4=-SfR/dxp;
      
    } else if(i==Ncells-1){	// Rightmost vertex
      TRACE(-1,"Rightmost vertex");

      d xim1=tube.geom.vx(i-1);
      d xi=tube.geom.vx(i);
      d dxm=xi-xim1;

      Whim1=-wLl;
      Whi=-wLr;
      Whip1=0;

      Wjim1=wLl;
      Wji=wLr;
      Wjip1=0;

      Wc1=-SfL/dxm;
      Wc2=SfL/dxm;
      Wc3=0;
      Wc4=0;

    } else{			// Normal interior vertex

      d xip1=tube.geom.vx(i+1);
      d xim1=tube.geom.vx(i-1);
      d xi=tube.geom.vx(i);
      d dxp=xip1-xi;
      d dxm=xi-xim1;
      
      Whim1=-wLl;
      Whi=wRl-wLr;
      Whip1=wRr;

      Wjim1=wLl;
      Wji=wLr-wRl;
      Wjip1=-wRr;
      
      Wc1=-SfL/dxm;
      Wc2=SfL/dxm;
      Wc3=SfR/dxp;
      Wc4=-SfR/dxp;
    }
    TRACE(0,"Energy constructor done");
  }
  vd Energy::Error(){		// Error in momentum equation
    TRACE(0,"Energy::Error()");
    vd error(Ns,fillwith::zeros);
    d gamma=this->gamma();

    vd Uti=vertex.U.tdata();
    vd pti=vertex.p.tdata();
    vd Ti=vertex.T.tdata();
    error+=vVf*DDTfd*vertex.p()/(gamma-1.0);
    error+=Whi*gamma*fDFT*(pti%Uti)/(gamma-1.0);
    error+=Wji*fDFT*(pti%Uti);
    error+=fDFT*(Wc2*kappaL()%Ti+Wc3*kappaR()%Ti);
    
    if(i>0){
      vd Utim1=left->U.tdata();
      vd ptim1=left->p.tdata();
      vd Tim1=left->T.tdata();
      error+=Whim1*gamma*fDFT*(ptim1%Utim1)/(gamma-1.0);
      error+=Wjim1*fDFT*(ptim1%Utim1);
      error+=Wc1*fDFT*(kappaL()%Tim1);
    }
    if(i<Ncells-1){
      vd Utip1=right->U.tdata();
      vd ptip1=right->p.tdata();
      vd Tip1=right->T.tdata();
      error+=Whip1*gamma*fDFT*(ptip1%Utip1)/(gamma-1.0);
      error+=Wjip1*fDFT*(ptip1%Utip1);
      error+=Wc4*fDFT*(kappaR()%Tip1);
    }

    // (Boundary source term)
    TRACE(-1,"vertex.esource called from Vertex object..."<<vertex.esource());
    
    error+=vertex.esource();

    return error;
  }
  dmat Energy::dpi(){
    TRACE(0,"Energy::dpi()");
    d T0=vertex.T(0);
    d gamma=tube.gas.gamma(T0);
    dmat dpi=zero;
    dpi+=(vVf/(gamma-1.0))*DDTfd;
    dpi+=Whi*(gamma/(gamma-1.0))*fDFT*diagtmat(vertex.U)*iDFT;
    return dpi;
  }

  dmat Energy::dUi(){
    TRACE(0,"Energy::dUi()");
    dmat dUi=zero;			    // Initialize with zeros
    d T0=vertex.T(0);
    d gamma=tube.gas.gamma(T0);
    dUi+=Whi*(gamma/(gamma-1.0))*fDFT*diagtmat(vertex.p)*iDFT;
    
    dUi+=Wji*fDFT*diagtmat(vertex.p)*iDFT;
    if(i>0)
      dUi+=Wjip1*fDFT*diagmat(left->p.tdata())*iDFT;
    if(i<Ncells-1)
      dUi+=Wjim1*fDFT*diagmat(right->p.tdata())*iDFT;
    return dUi;
  }
  dmat Energy::dTi(){
    TRACE(0,"Energy::dTi()");
    dmat dTi=zero;
    dTi+=fDFT*(Wc2*diagmat(kappaL())+Wc3*diagmat(kappaR()))*iDFT;
    return dTi;
  }
  dmat Energy::dpip1(){
    TRACE(0,"Energy::dpip1()");
    vd Uti=vertex.U.tdata();
    dmat dpip1=zero;
    d gamma=this->gamma();

    dpip1+=Wjip1*fDFT*diagmat(Uti)*iDFT;
    if(i<Ncells-1){
      vd Utip1=right->U.tdata();
      dpip1+=Whip1*(gamma/(gamma-1))*fDFT*diagmat(Utip1)*iDFT;
    }
    return dpip1;
  }
  dmat Energy::dUip1(){
    TRACE(0,"Energy::dUip1()");
    d gamma=this->gamma();
    dmat dUip1=zero;
    if(i<Ncells-1){
      vd ptip1=right->p.tdata();
      dUip1+=Whip1*(gamma/(gamma-1.0)*fDFT*diagmat(ptip1))*iDFT;
    }
    return dUip1;
  }
  dmat Energy::dTip1(){
    TRACE(0,"Energy::dTip1()");
    dmat dTip1=zero;
    if(i<Ncells-1)    
      dTip1+=Wc4*fDFT*diagmat(kappaR())*iDFT;
    return dTip1;
  }
  dmat Energy::dpim1(){
    TRACE(0,"Energy::dpim1()");
    d gamma=this->gamma();

    vd Uti=vertex.U.tdata();
    dmat dpim1=zero;
    if(i>0){
      vd Utim1=left->U.tdata();
      dpim1+=Whim1*(gamma/(gamma-1))*fDFT*diagmat(Utim1)*iDFT;
      dpim1+=Wjim1*fDFT*diagmat(Uti)*iDFT;
    }
    return dpim1;
  }
  dmat Energy::dUim1(){
    TRACE(0,"Energy::dUim1()");
    d gamma=this->gamma();

    dmat dUim1=zero;
    if(i>0){
      vd ptim1=left->p.tdata();
      dUim1+=Whim1*(gamma/(gamma-1.0)*fDFT*diagmat(ptim1))*iDFT;
    }
    return dUim1;
    }
  dmat Energy::dTim1(){
    TRACE(0,"Energy::dTim1()");
    dmat dTim1=zero;
    if(i>0)    
      dTim1+=Wc1*fDFT*diagmat(kappaL())*iDFT;
    return dTim1;
  }    
  vd Energy::kappaL(){
    TRACE(0,"Energy::kappaL()");
    
    vd kappaL(Ns,fillwith::zeros);
    vd Tti=vertex.T.tdata();
    vd kappait=tube.gas.kappa(Tti);

    if(i==0){
      vd Ttip1=right->T.tdata();
      vd kappaitp1=tube.gas.kappa(Ttip1);
      kappaL=wL1*kappaitp1+wL0*kappait;
    }
    else{
      vd Ttim1=left->T.tdata();
      vd kappaitm1=tube.gas.kappa(Ttim1);
      kappaL=wLr*kappait+wLl*kappaitm1;
      }

    return kappaL;
  }
  vd Energy::kappaR(){		// Returns thermal conductivity time domain data
    TRACE(0,"Energy::kappaR()");
    
    vd kappaR(Ns,fillwith::zeros);
    vd Tti=vertex.T.tdata();
    vd kappait=tube.gas.kappa(Tti);

    if(i==Ncells-1){
      vd Ttim1=left->T.tdata();
      vd kappaitm1=tube.gas.kappa(Ttim1);
      kappaR=wRNm2*kappaitm1+wRNm1*kappait;
    }
    else{
      vd Ttip1=right->T.tdata();
      vd kappaitp1=tube.gas.kappa(Ttip1);
      kappaR=wRl*kappait+wRr*kappaitp1;
      }
    return kappaR;
  }
  d Energy::gamma(){
    d T0=vertex.T(0);
    return tube.gas.gamma(T0);
  }
  Energy::~Energy(){
    TRACE(-5,"Energy destructor");
  }
    
} // namespace tube




