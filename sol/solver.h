#pragma once
#ifndef _SOLVER_H_
#define _SOLVER_H_
#include "system.h"

namespace tasystem{

  class Solver
  {
  public:

  public:
    TAsystem* sys=NULL;
    
    Solver(const TAsystem& tasys);
    Solver(const Solver& other);
    Solver& operator=(const Solver& other);

    void DoIter(d dampfac=1.0);
    ~Solver();
    
  };

} // namespace tasystem

#endif /* _SOLVER_H_ */
