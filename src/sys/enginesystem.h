#pragma once
#ifndef _ENGINESYSTEM_H_
#define _ENGINESYSTEM_H_

#include "pickadof.h"
#include "tasystem.h"

namespace tasystem{

  // The EngineSystem uses the given Globalconf omg as its guess
  // value. It will iteratively change the base frequency to solve the
  // system.
  class TripletList;
  
  class EngineSystem:public TaSystem{
    d mass_=0;
    PickADof tc;		    // Timing constraint
    PickADof av;		    // Reference amplitude value
    d getInitialMassFromGc() const; //
  private:
    TripletList Ljac(d dampfac);
    TripletList Mjac(d dampfac);    
    evd errorL();
    evd errorM();    
  public:    
    EngineSystem(const Globalconf& gc);
    EngineSystem(const EngineSystem&);
    EngineSystem(const TaSystem&); // Copy from base class
    EngineSystem& operator=(const EngineSystem&);
    virtual TaSystem* copy() const {return new EngineSystem(*this);}

    void setTimingConstraint(us segnr,us vertexnr,us varnr,us freqnr);
    void setAmplitudeDof(us segnr,us vertexnr,us varnr,us freqnr);    
    // Solving methods
    virtual esdmat jac(d dampfac=-1);
    virtual evd error();
    virtual void setRes(const vd& res);
    virtual evd getRes();

    // Later on, these two should be moved to private
    vd dmtotdx() const;
    vd domg();
    
    void setMass(d mass){mass_=mass;}
    d getMass() const {return mass_;}


    virtual void init();
    virtual void show(us detailnr=0);
    // If user decides to set a "custom mass" in the system, which does
    // not correspond to rho0*TotalFluidVolume
  };

} // namespace tasystem


#endif /* _ENGINESYSTEM_H_ */


