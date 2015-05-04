#pragma once
#ifndef _ENGINESYSTEM_H_
#define _ENGINESYSTEM_H_

#include "tasystem.h"

namespace tasystem{

  #ifdef SWIG
  %catches(std::exception,...) EngineSystem::EngineSystem(const Globalconf& gc);
  %catches(std::exception,...) EngineSystem::EngineSystem(const EngineSystem&);
  %catches(std::exception,...) EngineSystem::init();
  %catches(std::exception,...) EngineSystem::show(us detailnr=0);
  #endif // SWIG

  // The EngineSystem uses the given Globalconf omg as its guess
  // value. It will iteratively change the base frequency to solve the
  // system.
  class TripletList;
  
  class EngineSystem:public TaSystem{
    // The Dof which should be used for the phase constraint
    int PhaseConstraintDof=-1;
    // Pointer to the segment which provide a phaseDof
    const segment::Seg* phaseDofSeg=NULL;

  public:    
    EngineSystem(const Globalconf& gc);
    EngineSystem(const EngineSystem&);
    EngineSystem(const TaSystem&); // Copy from base class
    EngineSystem& operator=(const EngineSystem&)=delete;
    virtual TaSystem* copy() const {return new EngineSystem(*this);}

    // Overloaded virtuals:
    vd Error();
    void setRes(const vd& res);
    vd getRes();

    void init();
    void show(us detailnr=0);
  private:
    // Overloaded virtual:
    TripletList jacTriplets(d dampfac);
    // Derivative of equations to frequency
    vd domg();
  };

} // namespace tasystem


#endif /* _ENGINESYSTEM_H_ */


