#include "tubeequation.h"
#include "tube.h"

namespace tube{

  inline double min(d x,d y)  {
    return x<=y? x : y;
  }
  inline double max(d x,d y)  {
    return x<=y? y : x;
  }
  
  TubeEquation::TubeEquation(TubeVertex& tgp):Equation(tgp),v(tgp)
  {
    TRACE(0,"TubeEquation constructor");

  }
  TubeEquation::TubeEquation(const TubeEquation& other):TubeEquation(other.v){
    TRACE(0,"TubeEquation copy constructor");
  }
  void TubeEquation::Init(const Globalconf& gc){
    TRACE(6,"TubeEquation::Init(gc)");
    Equation::Init(gc);
    zero=zeros<dmat>(gc.Ns,gc.Ns);
  }
  dmat TubeEquation::diagtmat(const variable::var& v){
    dmat result(v.gc->Ns,vertex.gc->Ns,fillwith::zeros);
    result.diag()=v.tdata();
    return result;
  }
  // vd TubeEquation::getp0(){
  //   TRACE(0,"TubeEquation::getp0()");
  //   vd p0(v.gc->Ns,fillwith::zeros);
  //   p0(0)=v.gc->p0;
  //   return p0;
  // }
  vd TubeEquation::getp0t(){
    TRACE(0,"TubeEquation::getp0t()");
    vd p0(v.gc->Ns,fillwith::ones);
    p0*=v.gc->p0;
    return p0;
  }
  dmat  TubeEquation::Jac(){
    // Compute the Jacobian for the subsystem around the current gridpoint

    TRACE(0,"TubeEquation::Jac()");
    const us Ns=v.gc->Ns;
    TRACE(2,"Assignment of Ns survived, Ns="<< v.gc->Ns);
    us bw=Ns-1;

    // For an unconnected boundary node, we need to shift all
    // equations one block to make space for connection to i-2 and i+2
    // derivatives
    dmat result(Ns,3*Neq*Ns,fillwith::zeros);
    // Order is: rho,U,T,p,Tw
    TRACE(-1,"Ns:" << Ns);
    // TRACE(-1,"rhoim1 size:"<< rhoim1);
    TRACE(-2,"gc dft size:"<< v.gc->fDFT.size());
    // submat: first row,first col,last row, last col
    long int offset=0;
    if(v.i==v.Ncells-1 && v.right==NULL)
      offset=Ns*Neq;
    if(v.i==0 && v.left==NULL){ // Most left node
      // TRACE(100,"First vertex not coupled to other left vertex");
      offset=-Ns*Neq;
      result.submat(0,10*Ns,bw,10*Ns+bw)=drhoip2();
      result.submat(0,11*Ns,bw,11*Ns+bw)=dUip2();
      result.submat(0,12*Ns,bw,12*Ns+bw)=dTip2();
      result.submat(0,13*Ns,bw,13*Ns+bw)=dpip2();
      result.submat(0,14*Ns,bw,14*Ns+bw)=dTsip2();
    }
    else{
      // TRACE(100,"First vertex IS coupled to other vertex");
      result.submat(0,offset     ,bw,offset+     bw)=drhoim1();
      result.submat(0,offset+1*Ns,bw,offset+  Ns+bw)=dUim1();
      result.submat(0,offset+2*Ns,bw,offset+2*Ns+bw)=dTim1();
      result.submat(0,offset+3*Ns,bw,offset+3*Ns+bw)=dpim1();
      result.submat(0,offset+4*Ns,bw,offset+4*Ns+bw)=dTsim1();
    }
    if(v.i==v.Ncells-1 && v.right==NULL){
      // TRACE(100,"Last vertex not coupled to other right vertex");
      result.submat(0,0   ,bw,     bw)=drhoim2();
      result.submat(0,1*Ns,bw,  Ns+bw)=dUim2();
      result.submat(0,2*Ns,bw,2*Ns+bw)=dTim2();
      result.submat(0,3*Ns,bw,3*Ns+bw)=dpim2();
      result.submat(0,4*Ns,bw,4*Ns+bw)=dTsim2();
    }
    else{
      // TRACE(100,"Last vertex IS coupled to other right vertex");
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
    // TRACE(2,"TubeEquation Jacobian created. Returning to vertex Jacobian.");
    // if(v.i==v.Ncells-1)
      // cout << "Last vertex Jacobian: for equation\n" << result;
    // if(v.i==0)
      // cout << "First vertex Jacobian: for equation\n" << result;

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
  vd eps(const vd& nu,d kappa){
    vd eps(nu.size());
    for(us i=0;i<eps.size(); i++)
      {
	eps(i)=max(min(0.5,kappa*nu(i)),kappa*1e-2);
      }
    return eps;
  }

  dmat TubeEquation::d_r(){
    const us Ns=v.gc->Ns;
    if(v.i==v.Ncells-1)
      return d_l();
    else {
      dmat Dr(Ns,Ns,fillwith::zeros);
      d rj=v.gc->c0;
      vd eps1=eps(nu(),v.gc->kappa);

      Dr.diag()=rj*eps1;
      return Dr;
    }
  }

  
  dmat TubeEquation::d_l(){
    const us Ns=v.gc->Ns;
    if(v.i==0)
      return d_r();
    else{
      dmat Dl(Ns,Ns,fillwith::zeros);
      d rj=v.gc->c0;
      vd eps1=eps(nu(),v.gc->kappa);

      Dl.diag()=rj*eps1;
      return Dl;
    }
  }
  vd TubeEquation::nu(){
    const d& Ns=v.gc->Ns;
    vd pi(Ns);
    vd pip1(Ns);
    vd pim1(Ns);    
    if(v.i>0 && v.i<v.Ncells-1){
      pi=v.p();
      pip1=v.right->p();
      pim1=v.left->p();
    } else if(v.i==0){
      pi=v.p();
      pim1=v.right->p();
      pip1=v.right->right->p();
    }
    else{
      pi=v.left->p();
      pip1=v.left->p();
      pim1=v.left->left->p();
    } // Last node
      // return ones<vd>(Ns);
    vd half=vd(Ns); half.fill(0.5);
    d denominator=3.0*v.gc->p0;//abs(pim1+2*pi+pip1);
    vd numerator=abs(pim1-2*pi+pip1);
    vd num_over_denom(Ns);
    for(us k=0;k<Ns;k++){
      // if(abs(denominator(k))<=1e-10) // Safe from division by zero??
	// num_over_denom(k)=0.5;
      // else
      
      num_over_denom(k)=numerator(k)/denominator;///denominator(k);
      // num_over_denom(k)=0.5;      
    }
  return num_over_denom;
  }
  
  TubeEquation::~TubeEquation(){TRACE(-5,"~TubeEquation()");}
} // namespace tube
