#include "tubeequation.h"
#include "tube.h"
#include "tubevertex.h"

namespace tube{

  TubeEquation::TubeEquation(const Tube& tube,TubeVertex& tgp):Equation(tube.gc),i(tgp.i),tube(tube),vertex(tgp),left(vertex.left),right(vertex.right),fDFT(gc.fDFT),iDFT(gc.iDFT),DDTfd(gc.DDTfd),geom(tube.geom),Ncells(geom.Ncells)
  {
    TRACE(0,"TubeEquation constructor");
    // Geometrical parameters


    zero=zeros<dmat>(tube.gc.Ns,tube.gc.Ns);
  }
  TubeEquation::TubeEquation(const TubeEquation& other):TubeEquation(other.tube,other.vertex){
    TRACE(0,"TubeEquation copy constructor");
  }
  dmat TubeEquation::diagtmat(const variable::var& v){
    dmat result(Ns,Ns,fillwith::zeros);
    result.diag()=v.tdata();
    return result;
  }
  vd TubeEquation::getp0(){
    TRACE(0,"TubeEquation::getp0()");
    vd p0(Ns,fillwith::zeros);
    p0(0)=tube.gc.p0;
    return p0;
  }
  vd TubeEquation::getp0t(){
    TRACE(0,"TubeEquation::getp0t()");
    vd p0(Ns,fillwith::ones);
    p0*=tube.gc.p0;
    return p0;
  }
  dmat  TubeEquation::Jac(){
    // Compute the Jacobian for the subsystem around the current gridpoint
    TRACE(0,"TubeEquation::Jac()");
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

  dmat TubeEquation::drhoim2(){
    TRACE(0,"TubeEquation::drhoim2()");
    return zero;}
  dmat TubeEquation::dUim2(){
    TRACE(0,"TubeEquation::dUim2()");
    return zero;}
  dmat TubeEquation::dTim2(){
    TRACE(0,"TubeEquation::dTim2()");
    return zero;}
  dmat TubeEquation::dpim2(){
    TRACE(0,"TubeEquation::dpim2()");
    return zero;}
  dmat TubeEquation::dTsim2(){
    TRACE(0,"TubeEquation::dTsim2()");
    return zero;}
  
  
  dmat TubeEquation::drhoim1(){
    TRACE(0,"TubeEquation::drhoim1()");
    return zero;}
  dmat TubeEquation::dUim1(){
    TRACE(0,"TubeEquation::dUim1()");
    return zero;}
  dmat TubeEquation::dTim1(){
    TRACE(0,"TubeEquation::dTim1()");
    return zero;}
  dmat TubeEquation::dpim1(){
    TRACE(0,"TubeEquation::dpim1()");
    return zero;}
  dmat TubeEquation::dTsim1(){
    TRACE(0,"TubeEquation::dTsim1()");
    return zero;}

  dmat TubeEquation::drhoi(){
    TRACE(0,"TubeEquation::drhoi()");
    return zero;}
  dmat TubeEquation::dUi(){
    TRACE(0,"TubeEquation::dUi()");
    return zero;}
  dmat TubeEquation::dTi(){
    TRACE(0,"TubeEquation::dTi()");
    return zero;}
  dmat TubeEquation::dpi(){
    TRACE(0,"TubeEquation::dpi()");
    return zero;}
  dmat TubeEquation::dTsi(){
    TRACE(0,"TubeEquation::dTsi()");
    return zero;}
  dmat TubeEquation::drhoip1(){
    TRACE(0,"TubeEquation::drhoip1()");
    return zero;}
  dmat TubeEquation::dUip1(){
    TRACE(0,"TubeEquation::dUip1()");
    return zero;}
  dmat TubeEquation::dTip1(){
    TRACE(0,"TubeEquation::dTip1()");
    return zero;}
  dmat TubeEquation::dpip1(){
    TRACE(0,"TubeEquation::dpip1()");
    return zero;}
  dmat TubeEquation::dTsip1(){
    TRACE(0,"TubeEquation::dTsip1()");
    return zero;}

  dmat TubeEquation::drhoip2(){
    TRACE(0,"TubeEquation::drhoip2()");
    return zero;}
  dmat TubeEquation::dUip2(){
    TRACE(0,"TubeEquation::dUip2()");
    return zero;}
  dmat TubeEquation::dTip2(){
    TRACE(0,"TubeEquation::dTip2()");
    return zero;}
  dmat TubeEquation::dpip2(){
    TRACE(0,"TubeEquation::dpip2()");
    return zero;}
  dmat TubeEquation::dTsip2(){
    TRACE(0,"TubeEquation::dTsip2()");
    return zero;}

  // Artificial viscosity matrices
  dmat TubeEquation::D_r(){
    const us& Ns=gc.Ns;		// Number of samples
    if(i==Ncells-1)
      return D_l();
    else {
      dmat Dr(Ns,Ns,fillwith::zeros);
      // Wesselings book: rj+0.5=speed of sound
      Dr.diag()=gc.c0*gc.kappa*nu()*gc.dx/vertex.dxp;
      return Dr;
    }
  }
  dmat TubeEquation::D_l(){
    const us& Ns=gc.Ns;		// Number of samples
    if(i==0)
      return D_r();
    else{
      dmat Dl(Ns,Ns,fillwith::zeros);
      // Wesselings book: rj+0.5=speed of sound
      Dl.diag()=gc.c0*gc.kappa*nu()*gc.dx/vertex.dxm;
      return Dl;
    }
  }
  vd TubeEquation::nu(){
    vd pi(Ns);
    vd pip1(Ns);
    vd pim1(Ns);    
    if(i>0&&i<Ncells-1){
      pi=vertex.p();
      pip1=right->p();
      pim1=left->p();
    } else if(i==0){
      pi=vertex.p();
      pim1=right->p();
      pip1=right->right->p();
    }
    else{
      pi=left->p();
      pip1=left->p();
      pim1=left->left->p();
    } // Last node
      // return ones<vd>(Ns);
    vd half=vd(Ns); half.fill(0.5);
    vd denominator=abs(pim1+2*pi+pip1);
    vd numerator=(pim1-2*pi+pip1);
    vd num_over_denom(Ns);
    for(us k=0;k<Ns;k++){
      // if(abs(denominator(k))<=1e-10) // Safe from division by zero??
	// num_over_denom(k)=0.01;
      // else
	// num_over_denom(k)=numerator(k)/denominator(k);
      num_over_denom(k)=0.1;
    }
    return min(half,num_over_denom);    
  }
  
  TubeEquation::~TubeEquation(){}
} // namespace tube
