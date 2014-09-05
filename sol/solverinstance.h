#ifndef _SOLVERINSTANCE_H_
#define _SOLVERINSTANCE_H_

namespace tasystem{
  class Solver;
  class SolverConfiguration;
  class SolverInstance{
    Solver* sol;
    SolverConfiguration* sc;
  public:
    SolverInstance(){}
    SolverInstance(Solver& sol,SolverConfiguration& sc);
    void operator()();
    void Kill();
};




} // namespace tasystem
#endif /* _SOLVERINSTANCE_H_ */
