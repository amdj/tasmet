// impedancebc.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "impedancebc.h"

namespace tube {
  using namespace tasystem;

  ImpedanceBc::ImpedanceBc(us segnr,Pos position,const var& z,d T0):
    TubeBc(segnr,position),
    z(z),
    T0(T0)
  {
    TRACE(15,"ImpedanceBc::ImpedanceBc()");
  }
  ImpedanceBc::ImpedanceBc(const ImpedanceBc& other,const TaSystem& sys):
    TubeBc(other,sys),
    z(other.z),
    T0(other.T0)
  {
    TRACE(15,"ImpedanceBc::ImpedanceBc(copy and init)");
    
  }
  void ImpedanceBc::updateNf() {
    z.updateNf();
  }
  void ImpedanceBc::setEqNrs(us firsteqnr) {
    TRACE(15,"void ImpedanceBc::setEqNrs()");
    TubeBc::setEqNrs(firsteqnr);
  } 
} // namespace tube

//////////////////////////////////////////////////////////////////////
