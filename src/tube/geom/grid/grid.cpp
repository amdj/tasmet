#include "grid.h"
#include "boundarylayer.h"
#define MAXGP 50000

namespace tube{

  void testgp(us& gp){
    if(gp<4 || gp >MAXGP)
      {
        WARN( "WARNING: Given number of gridpoints is "	\
	     << gp<<", which is too small, or larger than MAXGP.\n"
              << "MAXGP is: " << MAXGP << "."
              "Number of gridpoints will be coerced to 4");
        gp=4;
      }
  }

 
  Grid::Grid(us gp,d L):gp(gp),L(L) {
    TRACE(15,"Grid::Grid()");
    testgp(this->gp);
    if(L<=0){
      WARN("Illegal length chosen. Will be coerced to unity.");
      this->L=1;
    }
      
    makex();
  }
  Grid::Grid(const Grid& o):
    gp(o.gp),L(o.L),x(o.x)
  {
    if(o.bL)
      this->bL=o.bL->copy();
    if(o.bR)
      this->bR=o.bR->copy();
  }
  Grid& Grid::operator=(const Grid& o){
    if(&o!=this){
      gp=o.gp;
      L=o.L;
      x=o.x;
      if(o.bL)
        setLeftBl(*o.bL);
      if(o.bR)
        setRightBl(*o.bR);
    }
    return *this;
  }
  Grid::~Grid(){
    delete bL;
    delete bR;
  }
  void Grid::setLeftBl(const BoundaryLayer& b){
    TRACE(15,"Grid::setLeftBl()");
    delete bL;
    bL=b.copy();
    makex();
  }
  void Grid::setRightBl(const BoundaryLayer& b){
    TRACE(15,"Grid::setRightBl()");
    delete bR;
    bR=b.copy();
    makex();
  }
  void Smooth(vd& x){
    TRACE(35,"Smooth()");
    WARN("Make this function better!");
    // Do something like take every odd point and put it in the middle
    us size=x.size();
    if(x.size()>1)
      for(us i=1;i<size-1;i++){
        x(i)=0.5*(x(i+1)+x(i-1));
        i++;
      }
  }
  void Grid::makex() {
    TRACE(15,"Grid::makex()");
    x=linspace(0,L,gp);
    if(bL)
      makeLeftBl();
    if(bR)
      makeRightBl();
    Smooth(x);
  }

  void truncate(vd& x,d L){
    TRACE(15,"truncate()");
    // Truncate a vector till last element is still smaller than L
      vd newx=x;
      us i=0;
      while(newx(i)<L && i<newx.size()-1)
        i++;
      if(i>1)
        x=newx.subvec(0,i-1);
      else
        x=vd(0,fillwith::zeros);
  }
  void Grid::makeLeftBl() {
    TRACE(15,"Grid::makeLeftBl()");
    // Obtain a boundary layer
    vd bl=bL->getx();
    truncate(bl,getL()/2);
    
    // TRACE(15,"bl: "<< bl);
    us nminorig=0;
    us xlast=x.size();
    {
      d bllast=bl(bl.size()-1);
      // TRACE(15,"bllast"<< bllast);
      while(x(nminorig)<=bllast && nminorig<(xlast-1)){
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
  void Grid::makeRightBl() {
    TRACE(15,"Grid::makeRightBl()");
    // Obtain a boundary layer
    vd bl=bR->getx();
    truncate(bl,getL()/2);
    us nmaxorig=x.size()-1;
    {
      d bllast=bl(bl.size()-1);
      while(L-x(nmaxorig)<bllast && nmaxorig>0){
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
    x=newx;
  }

} // namespace tube
