#pragma once
#ifndef _SOLVER_H_
#define _SOLVER_H_
#include "tasystem.h"
#include "vtypes.h"
#include <memory>

#define SOLVER_MAXITER 100
namespace tasystem{
  using std::tuple;
  namespace boost{
    class thread;
  }
  class Solver
  {
    std::unique_ptr<TaSystem> tasystem;
    boost::thread *SolverThread;
  public:
    TaSystem& sys() const {return *tasystem.get();}
    Solver(const TaSystem& tasys);
    Solver(const Solver& other);
    Solver& operator=(const Solver& other);
    void solve(us maxiter=100,d funer=1e-8,d reler=1e-6,d dampfac=1.0);
    tuple<d,d> doIter(d dampfac=1.0);
  };

} // namespace tasystem

#endif /* _SOLVER_H_ */
