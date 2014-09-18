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
#include "energyeq.h"
#include "isentropiceq.h"
#include "stateeq.h"

namespace tube{

  // TwImpendance boundary condition can only be used for a right
  // side.
  
  
  
  class RightTwImpedanceEq: public TubeEquation{
  public:
    virtual TubeEquation* copy() const {return new RightTwImpedanceEq(*this);}
    virtual vd error(const TubeVertex&) const;
    virtual JacRow jac(const TubeVertex& v) const;
    virtual JacCol dUi(const TubeVertex&) const;
    virtual JacCol dUim1(const TubeVertex&) const;
    virtual JacCol dpR(const TubeVertex&) const;
  }; 

  class TwImpedance:public TubeBcVertex // Adiabatic impedance boundary condition
  {
  public:
    variable::var pr;    
    RightTwImpedanceEq twright; // Completely adjusted equation
    State midstate;
    // TwImpedanceEnergyEq   eright; // Completely adjusted equation
    // Isentropic is;
    virtual const variable::var& pR() const {return pr;}
    TwImpedance();
    TwImpedance(const TwImpedance& o);
    TwImpedance& operator=(const TwImpedance&);
    ~TwImpedance(){}
    virtual vd esource() const;		// Source term for constant temperature
    virtual void initTubeVertex(us i,const Tube& thisseg);
    virtual string getType() const {return string("TwImpedance");}
    virtual enum connectpos connectPos() const {return connectpos::right;}
    virtual TubeBcVertex* copy() const {return new TwImpedance(*this);}
  private:
    virtual void updateW(const SegBase&);
  };

} // namespace tube


#endif /* _TWIMPEDANCE_H_ */



