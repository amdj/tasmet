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
    PickADof tc;		    // Timing constraint
    PickADof av;		    // Reference amplitude value
  private:
    TripletList Ljac(d dampfac);
    TripletList Mjac(d dampfac);    
    vd errorL();
    vd errorM();    
  public:    
    EngineSystem(const Globalconf& gc);
    EngineSystem(const EngineSystem&);
    EngineSystem(const TaSystem&); // Copy from base class
    EngineSystem& operator=(const EngineSystem&)=delete;
    virtual TaSystem* copy() const {return new EngineSystem(*this);}

    void setTimingConstraint(us segnr,us cellnr,us Varnr,us freqnr);
    void setAmplitudeDof(us segnr,us cellnr,us Varnr,us freqnr);    

    // Overloaded virtuals:
    arma::sp_mat jac(d dampfac=-1);
    vd Error();
    void setRes(const vd& res);
    vd getRes();
    // Later on, these two should be moved to private
    vd domg();
    

    virtual void init();
    virtual void show(us detailnr=0);
    // If user decides to set a "custom mass" in the system, which does
    // not correspond to rho0*TotalFluidVolume
  };

} // namespace tasystem


#endif /* _ENGINESYSTEM_H_ */


