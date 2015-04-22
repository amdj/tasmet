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
    bool wait=true;
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
      wait=sc.wait;
      return *this;
    }
    #endif
    void setWait(bool w){wait=w;}
    void setFuntol(d ft){funtol=ft;}
    void setReltol(d ft){reltol=ft;}
    void setDampfac(d df){dampfac=df;}
    void setMindampfac(d df){mindampfac=df;}    
    void setMaxdampfac(d df){maxdampfac=df;}    
    SolverConfiguration(const SolverConfiguration& o){
      (*this)=o;
    }
    SolverConfiguration(bool wait1):SolverConfiguration()
    { wait=wait1;}
    SolverConfiguration(us maxiter=100,d funtol=1e-6,d reltol=1e-6,d mindampfac=1,d maxdampfac=1,bool wait=true):
      wait(wait),
      maxiter(maxiter),
      funtol(funtol),
      reltol(reltol),
      dampfac(maxdampfac),
      mindampfac(mindampfac),
      maxdampfac(maxdampfac)
    {}
    void show() const {
      const char* dowait=wait?"yes":"no";
      cout << "Solverconfiguration initialized.\n"
           << "Max. iterations: " << maxiter << "\n"                     \
           << "Funtol: " << funtol << "\n"                              \
           << "Reltol: " << reltol << "\n"
           << "Mindampfac: " << mindampfac << "\n"
           << "Maxdampfac: " << maxdampfac << "\n"
           << "Waiting for solution: " << dowait << ".\n";
      
    }
    
  };
} // namespace tasystem

#endif /* _SOLVERCONFIGURATION_H_ */
