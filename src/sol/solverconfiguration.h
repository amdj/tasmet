#ifndef _SOLVERCONFIGURATION_H_
#define _SOLVERCONFIGURATION_H_

#include <atomic>
#include "vtypes.h"

namespace tasystem{
  #ifndef SWIG
  SPOILNAMESPACE
  #endif

  class SolverConfiguration {
  public:
    #ifndef SWIG
    std::atomic<us> maxiter;
    std::atomic<d> funtol;
    std::atomic<d> reltol;
    std::atomic<d> dampfac;
    std::atomic<d> mindampfac;
    std::atomic<d> maxdampfac;    

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
    #endif
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
    }
    
  };
} // namespace tasystem

#endif /* _SOLVERCONFIGURATION_H_ */
