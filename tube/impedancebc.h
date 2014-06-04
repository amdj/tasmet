// file: bcvertex.h, created March 20th, 2014
// Author: J.A. de Jong

// bcvertex.h: external boundary conditions for tubes. This file
// contains the implementation of typical external boundary conditions
// for tubes as a custom vertex. Examples are adiabatic walls, isothermal walls and an
// adiabatic open pressure boundary conditions.
#pragma once

#ifndef _IMPEDANCEBC_H_
#define _IMPEDANCEBC_H_

#include "bcvertex.h"
#include "momentumeq.h"


namespace tube{

  class RightImpedanceMomentumEq:public Momentum{
  public:
    RightImpedanceMomentumEq(const Tube&,TubeBcVertex&,vd& Z);
    ~RightImpedanceMomentumEq(){}
    vd Error();
    dmat dUi();
    dmat dUim1();
    dmat dpi();
    dmat dpim1();
    dmat drhoim1();
    dmat drhoi();
    vd& Z;
    
  }; 

  class RightImpedance:public TubeBcVertex // Adiabatic impedance boundary condition
  {
  public:
    RightImpedance(const Tube& t,vd Z);
    ~RightImpedance(){}
    void updateW();
    vd Z;			// The impedance
    RightImpedanceMomentumEq mright; // Completely adjusted equation

  };



} // namespace tube


#endif /* _IMPEDANCEBC_H_ */
