#pragma once
#ifndef _SOLVER_H_
#define _SOLVER_H_
#include "solverconfiguration.h"
#include "vtypes.h"
#include <memory>
#include <boost/thread.hpp>
#include <boost/atomic.hpp>
#define SOLVER_MAXITER 100
namespace tasystem{
  using std::tuple;
  class TaSystem;
  class Solver
  {
    std::unique_ptr<TaSystem> tasystem;
    std::unique_ptr<boost::thread> solverThread;
    std::unique_ptr<SolverConfiguration> sc;
    boost::atomic<d> dampfac;

  public:
    Solver(const TaSystem& tasys);
    TaSystem& sys() { return *tasystem.get();}
    Solver(const Solver& other);
    Solver& operator=(const Solver& other);
    void stop();
    void solve(us maxiter=0,d funer=1e-8,d reler=1e-6,d dampfac=1.0,bool wait=true);
    tuple<d,d> doIter(d dampfac=-1.0);
    ~Solver();
  };

} // namespace tasystem

#endif /* _SOLVER_H_ */
