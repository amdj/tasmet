// file: bcvertex.h, created March 20th, 2014
// Author: J.A. de Jong

// bcvertex.h: external boundary conditions for tubes. This file
// contains the implementation of typical external boundary conditions
// for tubes as a custom vertex. Examples are adiabatic walls, isothermal walls and an
// adiabatic open pressure boundary conditions.
#pragma once
#ifndef _TWIMPEDANCE_H_
#define _TWIMPEDANCE_H_

#include "tubebcvertex.h"
#include "momentumeq.h"


namespace tube{
  using segment::connectpos;

  class TwImpedanceMomentumEq:public Momentum{
  public:
    TwImpedanceMomentumEq(TubeBcVertex&);
    ~TwImpedanceMomentumEq(){}
    vd Error();
    dmat dUi();
    dmat dUim1();
    dmat dpi();
    dmat dpim1();
    dmat drhoim1();
    dmat drhoi();
  }; 

  class TwImpedance:public TubeBcVertex // Adiabatic impedance boundary condition
  {
  public:
    TwImpedanceMomentumEq mright; // Completely adjusted equation
    
    TwImpedance(us segnr);
    TwImpedance(const TwImpedance& o);
    TwImpedance& operator=(const TwImpedance&);
    ~TwImpedance(){}
    virtual vd esource();		// Source term for constant temperature
    virtual void Init(us i,const SegBase& thisseg);
    virtual string gettype() const {return string("TwImpedance");}
    virtual enum connectpos connectPos() const {return connectpos::right;}
  private:
    virtual void updateW(const SegBase&);
    d xhalf;			// Distance to right wall
    
  };



} // namespace tube


#endif /* _TWIMPEDANCE_H_ */



