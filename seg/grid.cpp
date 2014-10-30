#include "grid.h"
#include "fsolve.h"
namespace segment{


  class err_alpha{
  public:
    d L_ov_dxb;
    us n;
    err_alpha(d L_ov_dxb,us n):L_ov_dxb(L_ov_dxb),n(n){assert(n>1);}
    d operator()(const d& alpha){
      return  L_ov_dxb-(pow(alpha,double(n))-1.0)/(alpha-1.0);
    }
  };
  
  d findalpha(d L_ov_dxb,us n){
    // Compute the grid expansion factor
    err_alpha err(L_ov_dxb,n);
    math_common::Fsolverd solver;
    // solver.setVerbose(true);
    math_common::dfun fun=std::bind(&err_alpha::operator(),&err,_1);    
    return solver(fun,2);
  }
  vd boundaryLayer(d dxb,d L,us n){
    TRACE(15,"boundaryLayer()");
    // Create a boundary layer of thickness L, with smallest grid
    // spacing dxb and n number of gridpoints
    assert(n>1);
    assert(dxb<L);
    d alpha=findalpha(L/dxb,n);
    TRACE(10,"alpha:"<< alpha);

    vd bl(n);
    bl(0)=0;
    for(us i=1;i<n;i++)
      bl(i)=bl(i-1)+dxb*pow(alpha,i-1);
    return bl;
  }

  Grid::Grid(us gp,d L):gp(gp),L(L) {
    TRACE(15,"Grid::Grid()");
    assert(gp>3);
    assert(L>0);
  }

  void Grid::setLeftBl(d dxb,d percentage,us n){
    TRACE(15,"Grid::setLeftBl()");
    assert(percentage>0 && percentage<50);
    dxbL=dxb;
    percL=percentage;
    nL=n;
    blleft=true;
    
  }
  void Grid::setRightBl(d dxb,d percentage,us n){
    TRACE(15,"Grid::setLeftBl()");
    assert(percentage>0 && percentage<50);
    dxbR=dxb;
    percR=percentage;
    nR=n;
    blright=true;
  }  
  vd Grid::getx() const{
    vd x=linspace(0,L,gp);
    if(blleft)
      makeLeftBl(x);
    if(blright)
      makeRightBl(x);
    return x;
  }
  
  void Grid::makeLeftBl(vd& x) const {
    d dxorig=x(1)-x(0);
    if(dxorig<dxbL){
      WARN("Boundary layer spacing wider than grid. Doing nothing");
      return;
    }
    const d blLmax=0.01*percL*x(x.size()-1);
    const vd bl=boundaryLayer(dxbL,blLmax,nL);

    // TRACE(15,"bl: "<< bl);
    us nminorig=0;
    {
      d bllast=bl(bl.size()-1);
      // TRACE(15,"bllast"<< bllast);
      while(x(nminorig)<bllast){
        nminorig++;
      }
    }
    vd newx(x.size()-nminorig+bl.size());
    us blsize=bl.size();
    us i=0;
    for(i=0;i<blsize;i++)
      newx(i)=bl(i);
    newx.subvec(bl.size(),newx.size()-1)=x.subvec(nminorig,x.size()-1);
    
    x=newx;
  }
  
 
  void Grid::makeRightBl(vd& x) const {
    TRACE(10,"Grid::setRightBl()");

    const d L=x(x.size()-1);
    const d dxorig=L-x(x.size()-2);
    if(dxorig<dxbR){
      WARN("Boundary layer spacing wider than grid. Doing nothing");
      return;
    }

    const d blLmax=0.01*percR*L;
    // Obtain a boundary layer
    const vd bl=boundaryLayer(dxbR,blLmax,nR);

    us nmaxorig=x.size()-1;
    {
      d bllast=bl(bl.size()-1);

      while(L-x(nmaxorig)<bllast){
        TRACE(15,"nmaxorig:"<< nmaxorig);
        nmaxorig--;
      }
    }
    vd newx(nmaxorig+1+bl.size());
    us blsize=bl.size();

    newx.subvec(0,nmaxorig)=x.subvec(0,nmaxorig);
    us i=0;
    for(;i<bl.size();i++)
      newx(newx.size()-1-i)=L-bl(i);

    // newx.subvec(minnumberold,newx.size()-1)=x;
    x=newx;
  }

}
