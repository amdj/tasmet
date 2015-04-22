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
#include "vtypes.h"
#include <memory>
#include <thread>
#include <cassert>
#include "vtypes.h"


namespace tasystem{

  struct ErrorVals{
    d funer;
    d reler;
    ErrorVals(d fe,d re): funer(fe),reler(re){}
    #ifndef SWIG                // Swig does not know about
                                // initializer lists
    ErrorVals(std::initializer_list<d> il)
    {
      assert(il.size()==2);
      auto it=il.begin();
      funer=*it++;
      reler=*it;
    }
    #endif
    ~ErrorVals(){ TRACE(20,"~ErrorVals()"); }
  };

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
  public:
    Solver(const SolverConfiguration sc=SolverConfiguration()) {
      SolverConfiguration::operator=(sc);
    }
    Solver(const Solver& other)=delete;    // No assignments as well
    Solver& operator=(const Solver& other)=delete;

    // Stop all solver threads
    void stop();

    #ifdef SWIG
    %newobject solve;
    #endif
    // Start solving a system
    void solve(TaSystem&);
    #ifndef SWIG
    void solve(TaSystem&,const SolverConfiguration& sc);
    #endif

    ~Solver();
  };

} // namespace tasystem

#endif /* _SOLVER_H_ */

//////////////////////////////////////////////////////////////////////
