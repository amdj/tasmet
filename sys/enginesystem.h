#pragma once
#ifndef _ENGINESYSTEM_H_
#define _ENGINESYSTEM_H_
#include "timingconstraint.h"
#include "tasystem.h"

namespace tasystem{

  // The EngineSystem uses the given Globalconf omg as its guess
  // value. It will iteratively change the base frequency to solve the
  // system.
  
  class EngineSystem:public TaSystem{
    TimingConstraint tc;
    d getInitialMassFromGc() const; //
    us getNVertex() const;
  public:    
    EngineSystem(const Globalconf& gc,TimingConstraint tc);
    EngineSystem(const EngineSystem&);
    EngineSystem& operator=(const EngineSystem&);
    virtual TaSystem* copy() const {return new EngineSystem(*this);}
    virtual esdmat jac();
    virtual void setRes(const vd& res);
    virtual evd getRes();
    vd dmtotdx() const;
    vd domg();
    virtual evd error();
    virtual void init();
    virtual void show(bool showvertices);
    // If user decides to set a "custom mass" in the system, which does
    // not correspond to rho0*TotalFluidVolume
    void setInitialMass(d mass) {gc.setMass(mass);} // Check for
						    // nonzero mass
  };

} // namespace tasystem


#endif /* _ENGINESYSTEM_H_ */


