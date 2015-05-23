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
  var adiabaticTemp(const var& pres){
    TRACE(10,"adiabaticTemp()");
    const Globalconf* gc=&pres.gc();
    d T0=gc->T0();
    d gamma=gc->gas().gamma(T0);
    vd p0(Ns,fillwith::ones); p0*=gc->p0();
    vd Tbct=T0*pow((p0+pres.tdata())/p0,(gamma-1.0)/gamma);		// Adiabatic compression/expansion
    var res(pres.gc());
    res.settdata(Tbct);
    return res;
  } 

  
} // namespace tasystem

//////////////////////////////////////////////////////////////////////
