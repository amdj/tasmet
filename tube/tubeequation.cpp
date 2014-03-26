#include "tubeequation.h"
#include "tube.h"
#include "vertex.h"


namespace tube{

  Equation::Equation(const Tube& tube,TubeVertex& tgp):i(tgp.i),tube(tube),vertex(tgp),left(vertex.left),right(vertex.right),vop(tube.vop),fDFT(vop.fDFT),iDFT(vop.iDFT),DDTfd(vop.DDTfd),Ns(vop.Ns),geom(tube.geom),Ncells(geom.Ncells)
  {
// Sf=tA
    TRACE(0,"Equation constructor");
    // Geometrical parameters

    TRACE(-1,"i:"<< i);
    SfL=geom.Sf(i);
    SfR=geom.Sf(i+1);
    TRACE(-1,"SfL,SfR created");
    vSf=geom.vSf(i);
    vSs=geom.vSs(i);
    vVf=geom.vVf(i);
    vVs=geom.vVs(i);

    wLl=0; wLr=0; wRr=0; wRl=0;		// Initialize weight functions to zero
    xR=tube.geom.x(i+1);		// Position of right cell wall
    xL=tube.geom.x(i);		// Position of left cell wall
    const vd& vx=tube.geom.vx;
    const d& vxi=vx(i);
    const us& Ncells=tube.geom.Ncells;
    d vxip1=0;
    d vxim1=0;

    if(i>0){
      vxim1=vx(i-1);
      wLl=(vxi-xL)/(vxi-vxim1);
      wLr=(xL-vxim1)/(vxi-vxim1);
    }
    if(i<Ncells-1){
      TRACE(-1,"i:"<<i);
      TRACE(0,"xR:"<<xR);
      TRACE(0,"vxi"<<vxi)
      vxip1=vx(i+1);
      wRr=(xR-vxi)/(vxip1-vxi);
      wRl=(vxip1-xR)/(vxip1-vxi);
    }
    // special weight function part
    wL0=wL1=wRNm1=wRNm2=0;	// Put these weight functions to zero
    if(i==0){
      wL0=vxip1/(vxip1-vxi);
      wL1=-vxi/(vxip1-vxi);
    }
    if(i==Ncells-1){
      wRNm1=(vxim1-xR)/(vxim1-vxi);
      wRNm2=(xR-vxi)/(vxim1-vxi);
    }
    // end special weight function part
    TRACE(-1,"vertex i:"<<i);
    TRACE(-1,"vertex density:"<<vertex.rho());
    

    // TRACE(0,"Ns:"<<tube.gc.Ns);
    zero=zeros<dmat>(tube.gc.Ns,tube.gc.Ns);
  }
  Equation::Equation(const Equation& other):Equation(other.tube,other.vertex){
    TRACE(0,"Equation copy constructor");
  }
  dmat Equation::diagtmat(variable::var& v){
    dmat result(Ns,Ns,fillwith::zeros);
    result.diag()=v.tdata();
    return result;
  }
  dmat  Equation::Jac(){
    // Compute the Jacobian for the subsystem around the current gridpoint
    TRACE(0,"Equation::Jac()");
    const variable::varoperations& vop=tube.vop; // Reference to variable operations
    us Ns=vop.Ns;		// Number of samples
    us bw=Ns-1;
    dmat result(Ns,3*Neq*Ns,fillwith::zeros);
    // Order is: rho,U,T,p,Tw
    TRACE(-1,"Ns:" << Ns);
    // TRACE(-1,"rhoim1 size:"<< rhoim1);
    TRACE(-2,"vop dft size:"<< vop.fDFT.size());
    // submat: first row,first col,last row, last col
    TRACE(-1,"Filling drhoim1");
    result.submat(0,  0 ,bw,     bw)=drhoim1();
    TRACE(-1,"Filling dUim1");
    result.submat(0,Ns,bw,  Ns+bw)=dUim1();
    TRACE(-1,"Filling dTim1");
    result.submat(0,2*Ns,bw,2*Ns+bw)=dTim1();
    TRACE(-1,"Filling dpim1");
    result.submat(0,3*Ns,bw,3*Ns+bw)=dpim1();
    TRACE(-1,"Filling dTsim1");
    result.submat(0,4*Ns,bw,4*Ns+bw)=dTsim1();


    TRACE(-1,"Filling drhoi");
    result.submat(0,5*Ns,bw,5*Ns+bw)=drhoi();
    TRACE(-1,"Filling dUi");
    result.submat(0,6*Ns,bw,6*Ns+bw)=dUi();
    TRACE(-1,"Filling dTi");
    result.submat(0,7*Ns,bw,7*Ns+bw)=dTi();
    TRACE(-1,"Filling dpi");
    result.submat(0,8*Ns,bw,8*Ns+bw)=dpi();
    TRACE(-1,"Filling dTsi");
    result.submat(0,9*Ns,bw,9*Ns+bw)=dTsi();

    TRACE(-1,"Filling drhoip1");
    result.submat(0,10*Ns,bw,10*Ns+bw)=drhoip1();
    TRACE(-1,"Filling dUip1");
    result.submat(0,11*Ns,bw,11*Ns+bw)=dUip1();
    TRACE(-1,"Filling dTip1");
    result.submat(0,12*Ns,bw,12*Ns+bw)=dTip1();
    TRACE(-1,"Filling dpip1");
    result.submat(0,13*Ns,bw,13*Ns+bw)=dpip1();
    TRACE(-1,"Filling dTsip1");
    result.submat(0,14*Ns,bw,14*Ns+bw)=dTsip1();
    TRACE(-1,"Equation Jac done");
    return result;
  }
  dmat Equation::drhoim1(){
    TRACE(-1,"Equation::drhoim1()");
    return zero;}
  dmat Equation::dUim1(){
    TRACE(0,"Equation::dUim1()");
    return zero;}
  dmat Equation::dTim1(){
    TRACE(0,"Equation::dTim1()");
    return zero;}
  dmat Equation::dpim1(){
    TRACE(0,"Equation::dpim1()");
    return zero;}
  dmat Equation::dTsim1(){
    TRACE(0,"Equation::dTsim1()");
    return zero;}

  dmat Equation::drhoi(){
    TRACE(0,"Equation::drhoi()");
    return zero;}
  dmat Equation::dUi(){
    TRACE(0,"Equation::dUi()");
    return zero;}
  dmat Equation::dTi(){
    TRACE(0,"Equation::dTi()");
    return zero;}
  dmat Equation::dpi(){
    TRACE(0,"Equation::dpi()");
    return zero;}
  dmat Equation::dTsi(){
    TRACE(0,"Equation::dTsi()");
    return zero;}
  dmat Equation::drhoip1(){
    TRACE(0,"Equation::drhoip1()");
    return zero;}
  dmat Equation::dUip1(){
    TRACE(0,"Equation::dUip1()");
    return zero;}
  dmat Equation::dTip1(){
    TRACE(0,"Equation::dTip1()");
    return zero;}
  dmat Equation::dpip1(){
    TRACE(0,"Equation::dpip1()");
    return zero;}
  dmat Equation::dTsip1(){
    TRACE(0,"Equation::dTsip1()");
    return zero;}

  Equation::~Equation(){}
} // namespace tube
