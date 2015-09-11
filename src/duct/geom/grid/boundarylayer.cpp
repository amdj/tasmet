#include "exception.h"
#include "boundarylayer.h"
#include "fsolve.h"
#include "grid.h"
#include <cassert>

namespace duct{

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

  us findn(d L_ov_dxb,d alpha){

    d arg=L_ov_dxb*(alpha-1)+1;
    d n_decimal=log(arg)/log(alpha);
    return (us) int(ceil(n_decimal));
  }

  BoundaryLayer::BoundaryLayer(d dxb,d L,us n):dxb(dxb),n(n){
    alpha=findalpha(L/dxb,n);
    assert(dxb<L);
    VARTRACE(25,n);
    TRACE(25,"found n:"<<findn(L/dxb,alpha));
  }
  BoundaryLayer::BoundaryLayer(d dxb,d L,d alpha):
    dxb(dxb),alpha(alpha)
  {
    TRACE(15,"Grid::setRightBl(growthfac)");
    assert(dxb<L);
    n=findn(L/dxb,alpha);
  }
  vd BoundaryLayer::getx() const {
    TRACE(15,"BoundaryLayer::getx()");
    // Create a boundary layer of thickness L, with smallest grid
    // spacing dxb and n number of gridpoints
    VARTRACE(15,n);
    vd bl(n);
    bl(0)=0;
    for(us i=1;i<n;i++)
      bl(i)=bl(i-1)+dxb*pow(alpha,i-1);
    return bl;
  }

  AutoBoundaryLayer::AutoBoundaryLayer(d dxb,d alpha,const Grid& g) throw(std::exception)
  {
    TRACE(15,"AutoBoundaryLayer::AutoBoundaryLayer()");
    this->dxb=dxb;
    this->alpha=alpha;
    if(dxb>g.getL()/2 || dxb<=1e-15)
      throw MyError("Illegal minimal boundary layer thickness");
    d dxo=g.getL()/(g.getgp()-1);
    n=(us) int(ceil(log(dxo/dxb)/log(alpha)+1));
    VARTRACE(35,n);
  }

} // namespace duct


