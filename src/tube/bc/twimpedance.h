// file: bccell.h, created March 20th, 2014
// Author: J.A. de Jong

// bccell.h: external boundary conditions for tubes. This file
// contains the implementation of typical external boundary conditions
// for tubes as a custom cell. Examples are adiabatic walls, isothermal walls and an
// adiabatic open pressure boundary conditions.
#pragma once
#ifndef _TWIMPEDANCE_H_
#define _TWIMPEDANCE_H_

#include "tubebccell.h"
#include "momentumeq.h"
#include "energyeq.h"
#include "isentropiceq.h"
#include "stateeq.h"
namespace tube{


  class TwImpedanceMomentumEq:public Momentum{
  public:
    virtual vd error(const Cell&) const;
    virtual JacCol dUi(const Cell&) const;
    virtual JacCol dUim1(const Cell&) const;
    virtual JacCol dpR(const Cell& v) const;
    virtual Equation* copy() const {return new TwImpedanceMomentumEq(*this);}
  }; 
  // class TwImpedanceEnergyEq:public Energy{
  // private:
  //   const TwImpedance& impedancecell;
  // public:
  //   TwImpedanceEnergyEq(TwImpedance&);
  //   ~TwImpedanceEnergyEq(){}
  //   virtual vd error(const Cell&) const;
  //   virtual dmat dUi(const Cell&) const;
  //   virtual dmat dUim1(const Cell&) const;
  // }; 
  class RightTwImpedanceEq: public Equation{
  public:
    virtual Equation* copy() const {return new RightTwImpedanceEq(*this);}
    virtual vd error(const Cell&) const;
    virtual JacRow jac(const Cell& v) const;
    virtual JacCol dUi(const Cell&) const;
    // virtual JacCol dUim1(const Cell&) const;
    virtual JacCol dpR(const Cell&) const;
    virtual JacCol dpL(const Cell&) const;    
  }; 

  class TwImpedance:public BcCell // Adiabatic impedance boundary condition
  {
  public:
    // TwImpedanceMomentumEq mright; // Completely adjusted equation
    // TwImpedanceEnergyEq   eright; // Completely adjusted equation
    RightTwImpedanceEq righttwimp;
    // Isentropic is;
    variable::var pr;
    // State s;
    StateR sr;
    
    const variable::var& pR() const {return pr;}
    TwImpedance();
    TwImpedance(const TwImpedance& o);
    TwImpedance& operator=(const TwImpedance&);
    ~TwImpedance(){}
    virtual vd esource() const;		// Source term for constant temperature
    virtual void initCell(us i,const Tube& thisseg);
    virtual string getType() const {return string("TwImpedance");}
    virtual enum connectpos connectPos() const {return connectpos::right;}
    virtual BcCell* copy() const {return new TwImpedance(*this);}
  private:
    virtual void updateW(const SegBase&);
    // friend class TwImpedanceMomentumEq;
    // friend class TwImpedanceEnergyEq;
    
  };



} // namespace tube


#endif /* _TWIMPEDANCE_H_ */



