#include "grid.h"
namespace segment{

  Grid::Grid(us gp,d L) {
    TRACE(15,"Grid::Grid()");
    assert(gp>4);
    assert(L>0);
    x=linspace(0,L,gp);
  }
  void Grid::setLeftBl(d dxmin,d alpha){
    assert(alpha>1);
    if(blleft){
      WARN("Left boundary layer already set! Returning");
      return;
    }
    d dxorig=x(1)-x(0);
    if(dxorig<dxmin){
      WARN("Boundary layer spacing wider than grid. Doing nothing");
      return;
    }
    vd xold=x;


    us nmax=floor(1+log(dxorig/dxmin)/log(alpha));
    d newdx=0;
    vector<d> leftbl(nmax);
    leftbl[0]=0;
    us n=1;
    for(us n=1;n<nmax;n++){
      leftbl[n]=leftbl[n-1]+pow(alpha,n)*dxmin;
    }

    blleft=true;
  }
  
  void Grid::setRightBl(d dxmin,d alpha){
    if(blright){
      WARN("Left boundary layer already set! Returning");
      return;
    }

  }

}
