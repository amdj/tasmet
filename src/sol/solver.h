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
  using std::tuple;
  class TaSystem;

  #ifndef SWIG
  evd solvesys_eigen(const esdmat& K,const evd& f);
  #endif

  class Solver
  {
    TaSystem* tasystem=NULL;
    std::unique_ptr<boost::thread> solverThread;
  public:
    #ifndef SWIG
    SolverConfiguration sc;
    #endif
    Solver(const TaSystem& tasys);
    TaSystem& sys() { return *tasystem;}
    Solver(const Solver& other);
    Solver& operator=(const Solver& other);
    void stop();
    void solve(us maxiter=5000,d funtol=1e-8,d reltol=1e-6,d mindampfac=1e-2,d maxdampfac=1,bool wait=true);
    #ifndef SWIG
    tuple<d,d> doIter(d dampfac=-1);
    #endif
    ~Solver();
  };

} // namespace tasystem

#endif /* _SOLVER_H_ */
