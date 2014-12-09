#include "tubeequation.h"
#include "tubevertex.h"


namespace tube{
  inline double min(d x,d y)  {
    return x<=y? x : y;
  }
  inline double max(d x,d y)  {
    return x<=y? y : x;
  }

  // Artificial viscosity matrices
  // vd eps(const vd& nu,d kappa) {
  //   TRACE(5,"eps()");
  //   vd eps(nu.size());
  //   for(us i=0;i<eps.size(); i++)
  //     {
  //   // eps(i)=0.5;
  //   // eps(i)=min(0.5,kappa*nu(i));
  //   eps(i)=max(min(0.5,kappa*nu(i)),kappa*1e-3);
  //     }
  //   // TRACE(50,"Max eps:"<< max(eps));
  //   return eps;
  // }

  // dmat TubeEquation::d_r() const {
  //   TRACE(3,"TubeEquation::d_r()");
  //   const us Ns=v.gc->Ns();
  //   if(!v.right())
  //     return d_l(v);
  //   else {
  //     dmat Dr(Ns,Ns,fillwith::zeros);
  //     d rj=v.gc->c0();
  //     vd eps1=eps(nu(v),v.gc->kappa);
  //     Dr.diag()=rj*eps1;
  //     return Dr;
  //   }
  // }

  
  // dmat TubeEquation::d_l() const {
  //   TRACE(3,"TubeEquation::d_l()");
  //   const us Ns=v.gc->Ns();
  //   if(!v.left)
  //     return d_r(v);
  //   else{
  //     dmat Dl(Ns,Ns,fillwith::zeros);
  //     d rj=v.gc->c0();
  //     vd eps1=eps(nu(v),v.gc->kappa);
  //     Dl.diag()=rj*eps1;
  //     return Dl;
  //   }
  // }
  // vd TubeEquation::nu() const {
  //   TRACE(3,"TubeEquation::nu()");
  //   const d& Ns=v.gc->Ns();
  //   TRACE(3,"SFSG"<<Ns);    
  //   vd pi(Ns);
  //   vd pip1(Ns);
  //   vd pim1(Ns);    

  //   if(!v.left()! && !v.right()){
  //     pi=v.pL()();
  //     pip1=v.right->p();
  //     pim1=v.left->p();
  //   } else if(v.i==0){
  //     pi=v.p();
  //     pim1=v.right->p();
  //     pip1=v.right->right->p();
  //   }
  //   else if(v.i==v.nCells-1){
  //     pi=v.left->p();
  //     pip1=v.left->p();
  //     pim1=v.left->left->p();
  //   }
  //   else{
  //     WARN("Impossible! v.i falls out of range. Something wrong. Aborting...");
  //     abort();
  //   }
  //   // Last node
  //     // return ones<vd>(Ns);
  //   vd half=vd(Ns); half.fill(0.5);
  //   TRACE(3,"SFSG");
  //   d denominator=3.0*v.gc->p0;//abs(pim1+2*pi+pip1);
  //   vd numerator=abs(pim1-2*pi+pip1);
  //   vd num_over_denom(Ns);
  //   for(us k=0;k<Ns;k++){
  //     num_over_denom(k)=numerator(k)/denominator;///denominator(k);
  //   }
  // return num_over_denom;
  // }

} // namespace tube
