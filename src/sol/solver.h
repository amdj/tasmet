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
#define SOLVER_MAXITER 100
namespace tasystem{
  using std::tuple;
  class TaSystem;

  evd solvesys_eigen(const esdmat& K,const evd& f);

  class Solver
  {
    TaSystem* tasystem=NULL;
    std::unique_ptr<boost::thread> solverThread;
    std::unique_ptr<SolverConfiguration> sc;

  public:
    Solver(const TaSystem& tasys);
    TaSystem& sys() { return *tasystem;}
    Solver(const Solver& other);
    Solver& operator=(const Solver& other);
    void stop();
    void solve(us maxiter=0,d funer=1e-8,d reler=1e-6,d mindampfac=1e-2,d maxdampfac=1,bool wait=true);
    tuple<d,d> doIter(d dampfac=1.0);
    ~Solver();
  };

} // namespace tasystem

#endif /* _SOLVER_H_ */
