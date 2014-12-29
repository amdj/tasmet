#ifndef _SOLVERCONFIGURATION_H_
#define _SOLVERCONFIGURATION_H_
#define SOLVER_MAXITER 100

#include <boost/atomic.hpp>
#include "vtypes.h"
namespace tasystem{
  SPOILNAMESPACE
  
  class SolverConfiguration {
  public:
    boost::atomic<us> maxiter;
    boost::atomic<d> funtol;
    boost::atomic<d> reltol;
    boost::atomic<d> dampfac;
    boost::atomic<d> mindampfac;
    boost::atomic<d> maxdampfac;    
    SolverConfiguration& operator=(const SolverConfiguration& sc){
      d maxiter=sc.maxiter;
      this->maxiter=maxiter;
      d funtol=sc.funtol;
      this->funtol=funtol;
      d reltol=sc.reltol;
      this->reltol=reltol;
      d dampfac=sc.dampfac;
      this->dampfac=dampfac;
      d mindampfac=sc.mindampfac;
      this->mindampfac=mindampfac;
      d maxdampfac=sc.maxdampfac;
      this->maxdampfac=maxdampfac;
      return *this;
    }
    SolverConfiguration(us maxiter=5000,d funtol=1e-6,d reltol=1e-6,d mindampfac=1,d maxdampfac=1):
      maxiter(maxiter),
      funtol(funtol),
      reltol(reltol),
      dampfac(maxdampfac),
      mindampfac(mindampfac),
      maxdampfac(maxdampfac)
    {
      cout << "Solverconfiguration initialized.\n"
           << "Max. iterations: " << maxiter << "\n"                     \
           << "Funtol: " << funtol << "\n"                              \
           << "Reltol: " << reltol << "\n"
           << "Mindampfac: " << mindampfac << "\n"  \
           << "Maxdampfac: " << maxdampfac << "\n";
      if(maxiter==0)
        maxiter=SOLVER_MAXITER;
    }
    
  };
} // namespace tasystem

#endif /* _SOLVERCONFIGURATION_H_ */
