#ifndef _SOLVERCONFIGURATION_H_
#define _SOLVERCONFIGURATION_H_
#include "vtypes.h"

namespace tasystem{
  typedef unsigned int us;
  typedef double d;
  
  class SolverConfiguration
  {
  public:
    us maxiter;
    d funtol;
    d reltol;
    d dampfac;
    SolverConfiguration(us maxiter,d funtol,d reltol,d dampfac):
      maxiter(maxiter),funtol(funtol),reltol(reltol),dampfac(dampfac){
      cout << "Solver started. Max. iterations: " << maxiter << "\n"    \
           << "Funtol: " << funtol << "\n"                              \
           << "Reltol: " << reltol << "\n";
    }
    
  };
}

#endif /* _SOLVERCONFIGURATION_H_ */
