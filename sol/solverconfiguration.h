#ifndef _SOLVERCONFIGURATION_H_
#define _SOLVERCONFIGURATION_H_
#include "vtypes.h"
#include <boost/atomic.hpp>
namespace tasystem{
  SPOILNAMESPACE
  typedef unsigned int us;
  typedef double d;
  
  class SolverConfiguration
  {
  public:
    boost::atomic<us> maxiter;
    boost::atomic<d> funtol;
    boost::atomic<d> reltol;
    boost::atomic<d> dampfac;
    boost::atomic<d> mindampfac;
    boost::atomic<d> maxdampfac;    
    

    SolverConfiguration(us maxiter,d funtol,d reltol,d mindampfac,d maxdampfac):
      maxiter(maxiter),
      funtol(funtol),
      reltol(reltol),
      dampfac(maxdampfac),
      mindampfac(mindampfac),
      maxdampfac(maxdampfac)
    {
      cout << "Solver started. Max. iterations: " << maxiter << "\n"    \
           << "Funtol: " << funtol << "\n"                              \
           << "Reltol: " << reltol << "\n"
           << "Mindampfac: " << mindampfac << "\n"  \
           << "Maxdampfac: " << maxdampfac << "\n";
    }
    
  };
}

#endif /* _SOLVERCONFIGURATION_H_ */
