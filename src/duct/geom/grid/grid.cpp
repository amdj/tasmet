#include "grid.h"
#include "exception.h"
#include "constants.h"
#include "staticmsg.h"

common::StaticMsg<> msg;

namespace duct{

  LinearGrid::LinearGrid(us gp,d L):
    gp(gp),
    L(L)
  {
    if(gp<constants::mingp)
      throw MyError(msg("Number of gridpoints is lower than minimum. Minimum is: %d",constants::mingp));
    else if(gp>constants::maxgp)
      throw MyError(msg("Maximum number of gridpoints exceeded. Maximum is: %d",constants::maxgp));
    if(L<=0){
      throw MyError("Illegal length chosen.");
    }
  }
  BlGrid::BlGrid(d L,d dxb,d dxmid):
    L(L),dxb(dxb),dxmid(dxmid)
  {
    TRACE(15,"BlGrid::BlGrid()");
  }

  vd BlGrid::getx() const {
    TRACE(15,"BlGrid::getx()");

    d delta=2*acosh(sqrt(dxmid/dxb));

    d dxdxi0=0.5*L*(delta/tanh(delta/2))*(1-tanh(delta/2));
    d dxi=dxb/dxdxi0;
    us N=ceil(1/dxi);

    vd xi=linspace(0,1,N);
    
    return L*0.5*(1+tanh(delta*(xi-0.5))/tanh(delta/2));

    
  }

} // namespace duct
