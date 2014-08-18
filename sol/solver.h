#pragma once
#ifndef _SOLVER_H_
#define _SOLVER_H_
#include "system.h"
#include <math_common.h>
#include <memory>

#define SOLVER_MAXITER 100
namespace tasystem{
  using std::tuple;

  class Solver
  {
    std::unique_ptr<taSystem> tasystem;
  public:
    taSystem& sys() const {return *tasystem.get();}
    Solver(const taSystem& tasys);
    Solver(const Solver& other);
    Solver& operator=(const Solver& other);
    void solve(us maxiter=0,d funer=1e-8,d reler=1e-6);
    tuple<d,d> doIter(d dampfac=1.0);
  };

} // namespace tasystem

#endif /* _SOLVER_H_ */
