#include "grid.h"
#include "fsolve.h"
namespace segment{

  Grid::Grid(us gp,d L) {
    TRACE(15,"Grid::Grid()");
    assert(gp>3);
    assert(L>0);
    x=linspace(0,L,gp);
  }

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

    err_alpha err(L_ov_dxb,n);
    math_common::Fsolverd solver;
    solver.setVerbose(true);
    math_common::dfun fun=std::bind(&err_alpha::operator(),&err,_1);    
    return solver(fun,2);
  }





  
  vd boundaryLayer(d dxb,d L,us n){
    TRACE(15,"boundaryLayer()");
    assert(n>1);
    assert(dxb<L);
    d alpha=findalpha(L/dxb,n);
    TRACE(15,"alpha:"<< alpha);

    vd bl(n);
    bl(0)=0;
    for(int i=1;i<n;i++)
      bl(i)=bl(i-1)+dxb*pow(alpha,i-1);
    return bl;
  }


  
  void Grid::setLeftBl(d dxb,d percentage,us n){
    TRACE(15,"Grid::setLeftBl()");
    assert(percentage>0 && percentage<50);
    if(blleft){
      WARN("Left boundary layer already set! Returning");
      return;
    }
    d dxorig=x(1)-x(0);
    if(dxorig<dxb){
      WARN("Boundary layer spacing wider than grid. Doing nothing");
      return;
    }
    const vd oldx=x;
    assert(dxorig>dxb);
    const d blLmax=0.01*percentage*x(x.size()-1);

    const vd bl=boundaryLayer(dxb,blLmax,n);

    // TRACE(15,"bl: "<< bl);
    us nminorig=0;
    {
      d bllast=bl(bl.size()-1);
      TRACE(15,"bllast"<< bllast);
      TRACE(15,"oldx:"<<oldx);
      // us oldxlast=oldx.size()-1;
      while(oldx(nminorig)<bllast){
        TRACE(15,"nminorig:"<< nminorig);
        nminorig++;
      }
    }
    vd newx(oldx.size()-nminorig+bl.size());
    us blsize=bl.size();
    us i=0;
    for(i=0;i<blsize;i++)
      newx(i)=bl(i);
    newx.subvec(bl.size(),newx.size()-1)=oldx.subvec(nminorig,oldx.size()-1);
    
    x=newx;
    blleft=true;
  }
  
 
  void Grid::setRightBl(d dxb,d percentage,us n){
    TRACE(10,"Grid::setRightBl()");
    
    assert(percentage>0 && percentage<50);
    if(blright){
      WARN("Right boundary layer already set! Returning");
      return;
    }
    const d L=x(x.size()-1);
    const d dxorig=L-x(x.size()-2);
    if(dxorig<dxb){
      WARN("Boundary layer spacing wider than grid. Doing nothing");
      return;
    }

    const vd oldx=x;
    TRACE(15,"oldx:"<<oldx);
    const d blLmax=0.01*percentage*L;
    // Obtain a boundary layer
    const vd bl=boundaryLayer(dxb,blLmax,n);

    us nmaxorig=oldx.size()-1;
    {
      d bllast=bl(bl.size()-1);

      while(L-oldx(nmaxorig)<bllast){
        TRACE(15,"nmaxorig:"<< nmaxorig);
        nmaxorig--;
      }
    }
    vd newx(nmaxorig+1+bl.size());
    us blsize=bl.size();

    newx.subvec(0,nmaxorig)=oldx.subvec(0,nmaxorig);
    us i=0;
    for(;i<bl.size();i++)
      newx(newx.size()-1-i)=L-bl(i);
    cout << "newx:" << newx;

    
    // newx.subvec(minnumberold,newx.size()-1)=x;
    x=newx;
    blright=true;
  }

}
