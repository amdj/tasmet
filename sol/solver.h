#pragma once
#ifndef _SOLVER_H_
#define _SOLVER_H_
#include "system.h"
#include <math_common.h>
#define SOLVER_MAXITER 100
namespace tasystem{
  using std::tuple;

  class Solver
  {

  public:
    TAsystem sys;
    Solver(const TAsystem& tasys);
    Solver(const Solver& other);
    Solver& operator=(const Solver& other);
    void solve(us maxiter=0,d funer=1e-8,d reler=1e-6);
    tuple<d,d> doIter(d dampfac=1.0);
  };

} // namespace tasystem

#endif /* _SOLVER_H_ */
