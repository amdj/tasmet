// varutils.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "varutils.h"
#define Ns (gc->Ns())

namespace tasystem {
  
  var coldTemp(const var& pres){
    return var(pres.gc(),pres.gc().T0());
  }
  var adiabaticTemp(const var& pres,d T0){
    TRACE(10,"adiabaticTemp()");
    const Globalconf* gc=&pres.gc();
    
    if(T0<=0)
      T0=gc->T0();
    d gamma=gc->gas().gamma(T0);
    vd p0(Ns,fillwith::ones); p0*=gc->p0();

    // Adiabatic compression/expansion
    var res(pres.gc(),T0*pow((p0+pres.tdata())/p0,(gamma-1.0)/gamma),false);
    return res;
  } 

  
} // namespace tasystem

//////////////////////////////////////////////////////////////////////
