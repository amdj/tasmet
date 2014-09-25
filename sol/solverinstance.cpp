#include "solverinstance.h"
#include "tasystem.h"

#include "solver.h"
#include "solverconfiguration.h"

namespace tasystem{
  typedef tuple<d,d> dtuple;


  SolverInstance::SolverInstance(Solver& sol,SolverConfiguration& sc):
    sol(&sol),sc(&sc)
  {
    TRACE(10,"SolverInstance()");
    // Starting StartSolving Thread
  }
  void SolverInstance::Kill(){
    TRACE(18,"SolverInstance::Kill()");
    
  }

  void SolverInstance::operator()(){
    TRACE(15,"SolverInstance::operator()");
    assert(sol!=NULL);

    TaSystem& sys=sol->sys();
    evd oldres=sys.getRes();
    d funer=1.0;
    // For sure, we do at least one iteration
    d reler=1.0;
    us nloop=0;

    if(sc->maxiter==0)
      sc->maxiter=SOLVER_MAXITER;

    while((funer>sc->funtol || reler>sc->reltol) && nloop<sc->maxiter)
      {
        try{
          dtuple ers=sol->doIter(sc->dampfac);
          funer=std::get<0>(ers);
          if(!(funer>0)){
            WARN("Function error: "<< funer << " . Quiting solving procedure.");
            break;
          }

          reler=std::get<1>(ers);
          cout << green <<  "Iteration: "<<nloop<<" , function error: "<<funer<<" , relative error:" << reler<< ".\n" << def;
          nloop++;
        }
        catch(int Error){
          cout << "Solver failed. Resetting result data to old state.\n";
          vd oldres_arma(oldres.data(),oldres.size(),false,false);
          sys.setRes(oldres_arma);
          return;
        }
      }
    if(nloop==sc->maxiter)
      WARN("Solver reached maximum number of iterations! Results might not be reliable!");
  }
  
  
  
} // namespace tasystem
