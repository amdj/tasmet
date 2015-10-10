// solver.h
//
// Author: J.A. de Jong 
//
// Description:
// When solve() is called, a copy of itself is created on the heap,
// for which a thread will be run. After the thread is done, the
// updated (solved) TaSystem will overwrite the old one.
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include "solverconfiguration.h"
#include <memory>
#include <thread>
#include <cassert>
#include "vtypes.h"
#include "errorvals.h"
#include "solprogress.h"

namespace tasystem{

  #ifndef SWIG
  class TaSystem;
  #endif  // ifndef SWIG
  #ifdef SWIG
  %catches(std::exception,...) Solver::solve(TaSystem&);
  %catches(std::exception,...) Solver::solve(TaSystem&,const SolverConfiguration&,bool wait=true);
  %catches(std::exception,...) doIter(TaSystem*,SolverConfiguration*);
  #endif  // ifdef SWIG
  // To do only one iteration
  ErrorVals doIter(TaSystem* sys,SolverConfiguration* sc=NULL);
  
  class Solver:public SolverConfiguration  {
    // Solverthread
    std::unique_ptr<std::thread> solverThread;
    SolProgress sp;
  public:
    Solver(const SolverConfiguration sc=SolverConfiguration()) {
      SolverConfiguration::operator=(sc);
    }
    #ifndef SWIG
    Solver(const Solver& other)=delete;    // No assignments as well
    Solver& operator=(const Solver& other)=delete;
    #endif // ifndef SWIG
    void setSc(const SolverConfiguration& sc){ SolverConfiguration::operator=(sc); }
    // Stop all solver threads
    void stop();
    #ifndef SWIG
    friend void doSolve(Solver* s,TaSystem* thesys,bool isThread);
    #endif

    #ifdef SWIG
    %newobject solve;
    #endif
    // Start solving a system
    SolProgress solve(TaSystem&);
    SolProgress solve(TaSystem&,const SolverConfiguration& sc);
    SolProgress getSp() const {return sp;}
    ~Solver();
  };

} // namespace tasystem

#endif /* _SOLVER_H_ */

//////////////////////////////////////////////////////////////////////
