#include "tubeequation.h"
#include "tubevertex.h"
#include "tube.h"

namespace tube{
  inline double min(d x,d y)  {
    return x<=y? x : y;
  }
  inline double max(d x,d y)  {
    return x<=y? y : x;
  }
  vd TubeEquation::domg(const TubeVertex& v) const{
    return vd(v.gc->Ns,fillwith::zeros);
  } 
  dmat TubeEquation::jac(const TubeVertex& v) const {
    // Compute the Jacobian for the subsystem around the current gridpoint
    TRACE(0,"TubeEquation::Jac()");
    us i=v.i;
    us nCells=v.nCells;
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
    if(i==nCells-1 && v.right==NULL)
      offset=Ns*Neq;
    if(i==0 && v.left==NULL){ // Most left node
      // TRACE(100,"First vertex not coupled to other left vertex");
      offset=-Ns*Neq;
      result.submat(0,10*Ns,bw,10*Ns+bw)=drhoip2(v);
      result.submat(0,11*Ns,bw,11*Ns+bw)=dUip2(v);
      result.submat(0,12*Ns,bw,12*Ns+bw)=dTip2(v);
      result.submat(0,13*Ns,bw,13*Ns+bw)=dpip2(v);
      result.submat(0,14*Ns,bw,14*Ns+bw)=dTsip2(v);
    }
    else{
      // TRACE(100,"First vertex IS coupled to other vertex");
      result.submat(0,offset     ,bw,offset+     bw)=drhoim1(v);
      result.submat(0,offset+1*Ns,bw,offset+  Ns+bw)=dUim1(v);
      result.submat(0,offset+2*Ns,bw,offset+2*Ns+bw)=dTim1(v);
      result.submat(0,offset+3*Ns,bw,offset+3*Ns+bw)=dpim1(v);
      result.submat(0,offset+4*Ns,bw,offset+4*Ns+bw)=dTsim1(v);
    }
    if(i==nCells-1 && v.right==NULL){
      // TRACE(100,"Last vertex not coupled to other right vertex");
      result.submat(0,0   ,bw,     bw)=drhoim2(v);
      result.submat(0,1*Ns,bw,  Ns+bw)=dUim2(v);
      result.submat(0,2*Ns,bw,2*Ns+bw)=dTim2(v);
      result.submat(0,3*Ns,bw,3*Ns+bw)=dpim2(v);
      result.submat(0,4*Ns,bw,4*Ns+bw)=dTsim2(v);
    }
    else{
      // TRACE(100,"Last vertex IS coupled to other right vertex");
      result.submat(0,offset+10*Ns,bw,offset+10*Ns+bw)=drhoip1(v);
      result.submat(0,offset+11*Ns,bw,offset+11*Ns+bw)=dUip1(v);
      result.submat(0,offset+12*Ns,bw,offset+12*Ns+bw)=dTip1(v);
      result.submat(0,offset+13*Ns,bw,offset+13*Ns+bw)=dpip1(v);
      result.submat(0,offset+14*Ns,bw,offset+14*Ns+bw)=dTsip1(v);
    }

    // And filling middle part...
    result.submat(0,offset+5*Ns,bw,offset+5*Ns+bw)=drhoi(v);
    result.submat(0,offset+6*Ns,bw,offset+6*Ns+bw)=dUi(v);
    result.submat(0,offset+7*Ns,bw,offset+7*Ns+bw)=dTi(v);
    result.submat(0,offset+8*Ns,bw,offset+8*Ns+bw)=dpi(v);
    result.submat(0,offset+9*Ns,bw,offset+9*Ns+bw)=dTsi(v);
    return result;
  }



  // Artificial viscosity matrices
  vd eps(const vd& nu,d kappa) {
    TRACE(10,"Temporarily overwriting eps to 0.5");
    vd eps(nu.size());
    for(us i=0;i<eps.size(); i++)
      {
	// eps(i)=0.5;
	// eps(i)=min(0.5,kappa*nu(i));
	eps(i)=max(min(0.5,kappa*nu(i)),kappa*1e-2);
      }
    // TRACE(50,"Max eps:"<< max(eps));
    return eps;
  }

  dmat TubeEquation::d_r(const TubeVertex& v) const {
    TRACE(3,"TubeEquation::d_r()");
    const us Ns=v.gc->Ns;
    if(v.right==NULL)
      return d_l(v);
    else {
      dmat Dr(Ns,Ns,fillwith::zeros);
      d rj=v.gc->c0;
      vd eps1=eps(nu(v),v.gc->kappa);
      Dr.diag()=rj*eps1;
      return Dr;
    }
  }

  
  dmat TubeEquation::d_l(const TubeVertex& v) const {
    TRACE(3,"TubeEquation::d_l()");
    const us Ns=v.gc->Ns;
    if(v.left==NULL)
      return d_r(v);
    else{
      dmat Dl(Ns,Ns,fillwith::zeros);
      d rj=v.gc->c0;
      vd eps1=eps(nu(v),v.gc->kappa);
      Dl.diag()=rj*eps1;
      return Dl;
    }
  }
  vd TubeEquation::nu(const TubeVertex& v) const {
    TRACE(3,"TubeEquation::nu()");
    const d& Ns=v.gc->Ns;
    TRACE(3,"SFSG"<<Ns);    
    vd pi(Ns);
    vd pip1(Ns);
    vd pim1(Ns);    

    if(v.left!=NULL && v.right!=NULL){
      pi=v.p();
      pip1=v.right->p();
      pim1=v.left->p();
    } else if(v.i==0){
      pi=v.p();
      pim1=v.right->p();
      pip1=v.right->right->p();
    }
    else if(v.i==v.nCells-1){
      pi=v.left->p();
      pip1=v.left->p();
      pim1=v.left->left->p();
    }
    else{
      WARN("Impossible! v.i falls out of range. Something wrong. Aborting...");
      abort();
    }
    // Last node
      // return ones<vd>(Ns);
    vd half=vd(Ns); half.fill(0.5);
    TRACE(3,"SFSG");
    d denominator=3.0*v.gc->p0;//abs(pim1+2*pi+pip1);
    vd numerator=abs(pim1-2*pi+pip1);
    vd num_over_denom(Ns);
    for(us k=0;k<Ns;k++){
      num_over_denom(k)=numerator(k)/denominator;///denominator(k);
    }
  return num_over_denom;
  }
  // ############################################################
  // Below, there is nothing more than an outline of the empty stencil
  dmat TubeEquation::drhoim2(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::drhoim2()");
    return v.zero;}
  dmat TubeEquation::dUim2(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dUim2()");
    return v.zero;}
  dmat TubeEquation::dTim2(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dTim2()");
    return v.zero;}
  dmat TubeEquation::dpim2(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dpim2()");
    return v.zero;}
  dmat TubeEquation::dTsim2(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dTsim2()");
    return v.zero;}
  
  
  dmat TubeEquation::drhoim1(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::drhoim1()");
    return v.zero;}
  dmat TubeEquation::dUim1(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dUim1()");
    return v.zero;}
  dmat TubeEquation::dTim1(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dTim1()");
    return v.zero;}
  dmat TubeEquation::dpim1(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dpim1()");
    return v.zero;}
  dmat TubeEquation::dTsim1(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dTsim1()");
    return v.zero;}

  
  dmat TubeEquation::drhoi(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::drhoi()");
    return v.zero;}
  dmat TubeEquation::dUi(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dUi()");
    return v.zero;}
  dmat TubeEquation::dTi(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dTi()");
    return v.zero;}
  dmat TubeEquation::dpi(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dpi()");
    return v.zero;}
  dmat TubeEquation::dTsi(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dTsi()");
    return v.zero;}

  
  dmat TubeEquation::drhoip1(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::drhoip1()");
    return v.zero;}
  dmat TubeEquation::dUip1(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dUip1()");
    return v.zero;}
  dmat TubeEquation::dTip1(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dTip1()");
    return v.zero;}
  dmat TubeEquation::dpip1(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dpip1()");
    return v.zero;}
  dmat TubeEquation::dTsip1(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dTsip1()");
    return v.zero;}

  dmat TubeEquation::drhoip2(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::drhoip2()");
    return v.zero;}
  dmat TubeEquation::dUip2(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dUip2()");
    return v.zero;}
  dmat TubeEquation::dTip2(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dTip2()");
    return v.zero;}
  dmat TubeEquation::dpip2(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dpip2()");
    return v.zero;}
  dmat TubeEquation::dTsip2(const TubeVertex& v) const {
    TRACE(0,"TubeEquation::dTsip2()");
    return v.zero;}  

} // namespace tube
