#pragma once
#ifndef _ENERGYEQ_H_
#define _ENERGYEQ_H_
#include "ductequation.h"

namespace duct{
  SPOILNAMESPACE
  class HeatSource;
  
  class Energy:public Equation
  {
    const Duct* t=nullptr;
    d Wddt=0,Wddtkin=0;

  public:
    Energy(const Cell& v):Equation(v){}
    d Htot() const;             // Total enthalpy flow through this node
    virtual void init();
    virtual tasystem::JacRow jac() const;
    virtual enum EqType getType() const { return EqType::Ene;}
    virtual void show() const; 
    virtual vd error() const;			// Error in Energy equation at node i
    virtual void domg(vd&) const;

    d gamma() const;			// Time-avg ratio of specific heats

    static vd extrapolateHeatFlow(const Cell&); 
    static tasystem::JacRow dExtrapolateHeatFlow(const Cell&);

    // Total enthalpy flow through left cell wall
    static vd mHl(const Cell&);
    // Total enthalpy flow through righ cell wall
    static vd mHr(const Cell&);
    static vd mEkinl(const Cell&);
    static vd mEkinr(const Cell&);
    static tasystem::JacRow dmHl(const Cell&);
    static tasystem::JacRow dmHr(const Cell&);

    static vd QL(const Cell&);              // Heat conduction trough left wall
    static vd QR(const Cell&);              // Heat conduction trough right wall
    static tasystem::JacRow dQL(const Cell&);
    static tasystem::JacRow dQR(const Cell&);
  private:


    // tasystem::JacRow dQL() const;
    // tasystem::JacRow dQR() const;

  };				// Energy

  class BcCell;
  class ExtrapolateEnthalpyFlow{
  public:
    static vd extrapolateEnthalpyFlow(const BcCell&);
    static tasystem::JacRow dExtrapolateEnthalpyFlow(const BcCell&);
  };

}      // namespace duct
#endif /* _ENERGYEQ_H_ */
