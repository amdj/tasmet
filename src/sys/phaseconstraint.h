// phaseconstraint.h
//
// Author: J.A. de Jong 
//
// Description:
// Provides a class with information about putting a phase constraint
// on a specific dependent variable.
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef PHASECONSTRAINT_H
#define PHASECONSTRAINT_H
#include "constants.h"
#include "vtypes.h"
#include "exception.h"

namespace tasystem {
  #ifdef SWIG
  %catches(std::exception,...) PhaseConstraint::PhaseConstraint(Varnr var,us freqnr=2,
                                                               Pos pos=Pos::left);
  #endif
  struct PhaseConstraint {
    Varnr var=Varnr::none;
    us freqnr=2;
    Pos pos;
    PhaseConstraint(Varnr var,us freqnr=2,Pos pos=Pos::left):
      var(var),freqnr(freqnr),pos(pos){
      if(freqnr<1)
        throw MyError("Illegal freqnr given. If you would like to "
                      "put the imaginary part if the first harmonic to zero"
                      ", please set freqnr=2.");
    }
    ~PhaseConstraint(){}
  };

  
} // namespace tasystem


#endif // PHASECONSTRAINT_H
//////////////////////////////////////////////////////////////////////
