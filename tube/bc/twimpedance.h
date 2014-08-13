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
namespace tube{


  class TwImpedanceMomentumEq:public Momentum{
  public:
    virtual vd error(const TubeVertex&) const;
    virtual dmat dUi(const TubeVertex&) const;
    virtual dmat dUim1(const TubeVertex&) const;
  }; 
  // class TwImpedanceEnergyEq:public Energy{
  // private:
  //   const TwImpedance& impedancevertex;
  // public:
  //   TwImpedanceEnergyEq(TwImpedance&);
  //   ~TwImpedanceEnergyEq(){}
  //   virtual vd error(const TubeVertex&) const;
  //   virtual dmat dUi(const TubeVertex&) const;
  //   virtual dmat dUim1(const TubeVertex&) const;
  // }; 

  class TwImpedance:public TubeBcVertex // Adiabatic impedance boundary condition
  {
  public:
    TwImpedanceMomentumEq mright; // Completely adjusted equation
    // TwImpedanceEnergyEq   eright; // Completely adjusted equation
    // Isentropic is;
    TwImpedance(us segnr);
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
    friend class TwImpedanceMomentumEq;
    // friend class TwImpedanceEnergyEq;
    
  };



} // namespace tube


#endif /* _TWIMPEDANCE_H_ */



