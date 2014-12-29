#ifndef _SOLVERINSTANCE_H_
#define _SOLVERINSTANCE_H_

namespace tasystem{
  class Solver;

  class SolverInstance{
    Solver* sol;
  public:
    SolverInstance(Solver& sol);
    void operator()();
    void Kill();
  };



} // namespace tasystem
#endif /* _SOLVERINSTANCE_H_ */
