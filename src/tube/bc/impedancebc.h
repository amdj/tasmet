// file: bcvertex.h, created March 20th, 2014
// Author: J.A. de Jong

// bcvertex.h: external boundary conditions for tubes. This file
// contains the implementation of typical external boundary conditions
// for tubes as a custom vertex. Examples are adiabatic walls, isothermal walls and an
// adiabatic open pressure boundary conditions.
#pragma once

#ifndef _IMPEDANCEBC_H_
#define _IMPEDANCEBC_H_

#include "tubebcvertex.h"
#include "momentumeq.h"


namespace tube{


  class RightImpedanceMomentumEq:public Momentum{
  public:
    RightImpedanceMomentumEq(TubeBcVertex&,vd& Z);
    ~RightImpedanceMomentumEq(){}
    virtual TubeEquation* copy(){ return new RightImpedanceMomentumEq(*this);}
    vd error(const TubeVertex&) const;
    JacCol dUi(const TubeVertex&) const;
    JacCol dUim1(const TubeVertex&) const;
    // dmat dpi();
    // dmat dpim1();
    // dmat drhoim1();
    // dmat drhoi();
    vd& Z;
  }; 

  class RightImpedance:public TubeBcVertex // Adiabatic impedance boundary condition
  {
  public:
    vd Z;			// The impedance
    RightImpedanceMomentumEq mright; // Completely adjusted equation
    
    
    RightImpedance(vd Z);
    RightImpedance(const RightImpedance& o);
    RightImpedance& operator=(const RightImpedance&);
    ~RightImpedance(){}
    virtual void initTubeVertex(us i,const Tube&);
    virtual string getType() const {return string("RightImpedance");}
    virtual Pos connectPos() const {return Pos::right;}
    virtual TubeBcVertex* copy() const {return new RightImpedance(*this);}
  protected:
    void updateW(const Tube&);

  };



} // namespace tube


#endif /* _IMPEDANCEBC_H_ */



