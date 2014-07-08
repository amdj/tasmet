#pragma once
#ifndef _SOLVER_H_
#define _SOLVER_H_
#include "system.h"

namespace tasystem{

class Solver
{
public:
  Solver(TAsystem& tasys);
  void DoIter(d dampfac=1.0);
private:
  TAsystem sys;
};

} // namespace tasystem

#endif /* _SOLVER_H_ */
