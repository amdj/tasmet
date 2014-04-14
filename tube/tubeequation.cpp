#include "tubeequation.h"
#include "tube.h"
#include "vertex.h"


namespace tube{

  Equation::Equation(const Tube& tube,TubeVertex& tgp):i(tgp.i),tube(tube),vertex(tgp),left(vertex.left),right(vertex.right),gc(tube.gc),fDFT(gc.fDFT),iDFT(gc.iDFT),DDTfd(gc.DDTfd),Ns(gc.Ns),geom(tube.geom),Ncells(geom.Ncells)
  {
    TRACE(0,"Equation constructor");
    // Geometrical parameters

    TRACE(-1,"i:"<< i);
    // Left and right cross-sectional area
    SfL=geom.Sf(i);
    SfR=geom.Sf(i+1);
    // Geometric parameters
    vSf=geom.vSf(i);
    vSs=geom.vSs(i);
    vVf=geom.vVf(i);
    vVs=geom.vVs(i);

    // Initialize weight functions to zero
    wLl=0; wLr=0; wRr=0; wRl=0;
    
    xR=tube.geom.x(i+1);		// Position of right cell wall
    xL=tube.geom.x(i);			// Position of left cell wall
    
    const vd& vx=tube.geom.vx;
    const d& vxi=vx(i);
    d vxip1=0;
    d vxim1=0;
    // Initialize distances to next node to zero
    dxm=dxp=0;

    if(i>0){
      vxim1=vx(i-1);
      dxm=vxi-vxim1;
      // Left weight functions
      wLl=(vxi-xL)/(vxi-vxim1);
      wLr=(xL-vxim1)/(vxi-vxim1);
    }
    if(i<Ncells-1){
      vxip1=vx(i+1);
      dxp=vxip1-vxi;
      // Right weight functions
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
    // TRACE(-1,"vertex i:"<<i);
    // TRACE(-1,"vertex density:"<<vertex.rho());
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
  vd Equation::getp0(){
    TRACE(0,"Equation::getp0()");
    vd p0(Ns,fillwith::zeros);
    p0(0)=tube.gc.p0;
    return p0;
  }
  vd Equation::getp0t(){
    TRACE(0,"Equation::getp0t()");
    vd p0(Ns,fillwith::ones);
    p0*=tube.gc.p0;
    return p0;
  }
  dmat  Equation::Jac(){
    // Compute the Jacobian for the subsystem around the current gridpoint
    TRACE(0,"Equation::Jac()");
    const tasystem::Globalconf& gc=tube.gc; // Reference to variable
					    // operations
    const us& Ns=gc.Ns;		// Number of samples
    us bw=Ns-1;
    const us Ncells=tube.geom.Ncells;
    // For an unconnected boundary node, we need to shift all
    // equations one block to make space for connection to i-2 and i+2
    // derivatives
    dmat result(Ns,3*Neq*Ns,fillwith::zeros);
    // Order is: rho,U,T,p,Tw
    TRACE(-1,"Ns:" << Ns);
    // TRACE(-1,"rhoim1 size:"<< rhoim1);
    TRACE(-2,"gc dft size:"<< gc.fDFT.size());
    // submat: first row,first col,last row, last col

    long int offset=0;
    if(i==Ncells-1 && vertex.right==NULL)
      offset=Ns*Neq;
    if(i==0 && vertex.left==NULL){ // Most left node
      offset=-Ns*Neq;
      result.submat(0,10*Ns,bw,10*Ns+bw)=drhoip2();
      result.submat(0,11*Ns,bw,11*Ns+bw)=dUip2();
      result.submat(0,12*Ns,bw,12*Ns+bw)=dTip2();
      result.submat(0,13*Ns,bw,13*Ns+bw)=dpip2();
      result.submat(0,14*Ns,bw,14*Ns+bw)=dTsip2();
    }
    else{
      result.submat(0,offset     ,bw,offset+     bw)=drhoim1();
      result.submat(0,offset+1*Ns,bw,offset+  Ns+bw)=dUim1();
      result.submat(0,offset+2*Ns,bw,offset+2*Ns+bw)=dTim1();
      result.submat(0,offset+3*Ns,bw,offset+3*Ns+bw)=dpim1();
      result.submat(0,offset+4*Ns,bw,offset+4*Ns+bw)=dTsim1();
    }
    if(i==Ncells-1 && vertex.right==NULL){
      result.submat(0,0   ,bw,     bw)=drhoim2();
      result.submat(0,1*Ns,bw,  Ns+bw)=dUim2();
      result.submat(0,2*Ns,bw,2*Ns+bw)=dTim2();
      result.submat(0,3*Ns,bw,3*Ns+bw)=dpim2();
      result.submat(0,4*Ns,bw,4*Ns+bw)=dTsim2();
    }
    else{
      result.submat(0,offset+10*Ns,bw,offset+10*Ns+bw)=drhoip1();
      result.submat(0,offset+11*Ns,bw,offset+11*Ns+bw)=dUip1();
      result.submat(0,offset+12*Ns,bw,offset+12*Ns+bw)=dTip1();
      result.submat(0,offset+13*Ns,bw,offset+13*Ns+bw)=dpip1();
      result.submat(0,offset+14*Ns,bw,offset+14*Ns+bw)=dTsip1();
    }

    // And filling middle part...
    result.submat(0,offset+5*Ns,bw,offset+5*Ns+bw)=drhoi();
    result.submat(0,offset+6*Ns,bw,offset+6*Ns+bw)=dUi();
    result.submat(0,offset+7*Ns,bw,offset+7*Ns+bw)=dTi();
    result.submat(0,offset+8*Ns,bw,offset+8*Ns+bw)=dpi();
    result.submat(0,offset+9*Ns,bw,offset+9*Ns+bw)=dTsi();


    return result;
  }

  dmat Equation::drhoim2(){
    TRACE(0,"Equation::drhoim2()");
    return zero;}
  dmat Equation::dUim2(){
    TRACE(0,"Equation::dUim2()");
    return zero;}
  dmat Equation::dTim2(){
    TRACE(0,"Equation::dTim2()");
    return zero;}
  dmat Equation::dpim2(){
    TRACE(0,"Equation::dpim2()");
    return zero;}
  dmat Equation::dTsim2(){
    TRACE(0,"Equation::dTsim2()");
    return zero;}
  
  
  dmat Equation::drhoim1(){
    TRACE(0,"Equation::drhoim1()");
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

  dmat Equation::drhoip2(){
    TRACE(0,"Equation::drhoip2()");
    return zero;}
  dmat Equation::dUip2(){
    TRACE(0,"Equation::dUip2()");
    return zero;}
  dmat Equation::dTip2(){
    TRACE(0,"Equation::dTip2()");
    return zero;}
  dmat Equation::dpip2(){
    TRACE(0,"Equation::dpip2()");
    return zero;}
  dmat Equation::dTsip2(){
    TRACE(0,"Equation::dTsip2()");
    return zero;}

  // Artificial viscosity matrices
  dmat Equation::D_r(){
    const us& Ns=gc.Ns;		// Number of samples
    if(i==Ncells-1)
      return D_l();
    else {
      dmat Dr(Ns,Ns,fillwith::zeros);
      // Wesselings book: rj+0.5=speed of sound
      Dr.diag()=gc.c0*gc.kappa*nu()*gc.dx/dxp;
      return Dr;
    }
  }
  dmat Equation::D_l(){
    const us& Ns=gc.Ns;		// Number of samples
    if(i==0)
      return D_r();
    else{
      dmat Dl(Ns,Ns,fillwith::zeros);
      // Wesselings book: rj+0.5=speed of sound
      Dl.diag()=gc.c0*gc.kappa*nu()*gc.dx/dxm;
      return Dl;
    }
  }
  vd Equation::nu(){
    vd pi(Ns);
    vd pip1(Ns);
    vd pim1(Ns);    
    if(i>0&&i<Ncells-1){
      pi=vertex.p();
      pip1=right->p();
      pim1=left->p();
    } else if(i==0){
      pim1=vertex.p();
      pi=right->p();
      pip1=tube.vvertex[2]->p();
    }
    else{
      pi=left->p();
      pip1=vertex.p();
      pim1=tube.vvertex[i-2]->p();
    } // Last node
      // return ones<vd>(Ns);
    vd half=vd(Ns); half.fill(0.5);
    vd denominator=abs(pim1+2*pi+pip1);
    vd numerator=((pim1-pi)+(pip1-pi)*dxm/dxp);
    vd num_over_denom(Ns);
    for(us k=0;k<Ns;k++){
      // if(abs(denominator(k))<=1e-10) // Safe from division by zero??
	// num_over_denom(k)=0;
      // else
      // num_over_denom(k)=numerator(k)/denominator(k);
      num_over_denom(k)=0.1;
    }
    return min(half,num_over_denom);    
  }
  
  Equation::~Equation(){}
} // namespace tube






