#pragma once
#ifndef _SOLVER_H_
#define _SOLVER_H_
#include "solverconfiguration.h"
#include "vtypes.h"
#include "arma_eigen.h"
#include <memory>
#include <boost/thread.hpp>
#include <boost/atomic.hpp>
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

  evd solvesys_eigen(const esdmat& K,const evd& f);
  #endif  // ifndef SWIG
  #ifdef SWIG
  %catches(std::exception,...) Solver::Solver(const TaSystem& sys);
  %catches(std::exception,...) Solver::Solver(const Solver& other);
  %catches(std::exception,...) Solver::solve(us maxiter=5000,d funtol=1e-8,\
                                             d reltol=1e-6,d mindampfac=1e-2,\
                                             d maxdampfac=1,bool wait=true);

  #endif  // ifdef SWIG
  class Solver
  {
    TaSystem* tasystem=nullptr;
    std::unique_ptr<boost::thread> solverThread;
  public:
    #ifndef SWIG
    SolverConfiguration sc;
    #endif  // ifndef SWIG
    // The best way to initialize a solver is by using a TaSystem to
    // work on.
    Solver(const TaSystem& tasys);
    Solver(const Solver& other);

    // Return a reference to the TaSystem
    TaSystem& sys() const { return *tasystem;}
    Solver& operator=(const Solver& other) =delete;
    // Stop all solver threads
    void stop();

    // Start a solver thread. 
    void solve(us maxiter=5000,d funtol=1e-8,d reltol=1e-6,d mindampfac=1e-2,d maxdampfac=1,bool wait=true);
    ErrorVals doIter(d dampfac=-1);
    ~Solver();
  };

} // namespace tasystem

#endif /* _SOLVER_H_ */

