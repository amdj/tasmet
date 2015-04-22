#include "solver.h"
#include "tasystem.h"
#include "vtypes.h"
#include "exception.h"

namespace tasystem{
  using arma::sp_mat;

  // A solver always contains a valid system.
  void Solver::stop() {
    TRACE(15,"Solver::stop()");
    if(solverThread){
      d curmaxiter=maxiter;
      maxiter=0;
      cout << "Waiting for solver to finish...\n";
      solverThread->join();
      solverThread.reset();
      maxiter=curmaxiter;
    }
  }
  void doSolve(Solver* s,TaSystem* thesys) {
    TRACE(15,"doSolve()");
    SolverConfiguration& sc=*s;
    assert(thesys && s);
    d funer=1.0;
    // For sure, we do at least one iteration
    d reler=1.0;
    us nloop=0;
    vd oldres=thesys->getRes();
    while((funer>sc.funtol || reler>sc.reltol) && nloop<sc.maxiter)    {
      try{
        ErrorVals ers=doIter(thesys,&sc);
        funer=ers.funer;
        reler=ers.reler;
        if(funer<0){
          WARN("Function error: "<< funer << " . Quiting solving procedure.");
          return;
        }
        cout << green <<  "Iteration: "<<nloop<<" , function error: "<<funer<<" , relative error:" << reler<< "."<< def <<"\n";
        nloop++;
      }
      catch(...){
        cout << "Solver failed, probably due to a numerical problem. TaSystem not updated.\n";
        thesys->setRes(oldres);
        return;
      }
    } // while
    if(nloop==sc.maxiter)
      WARN("Solver reached maximum number of iterations! Results might not be reliable!");
    if(sc.maxiter==0)
      WARN("Solver stopped externally");
    cout << "Solver done.\n";
  }

  void Solver::solve(TaSystem& sys){
    TRACE(20,"Solver::solve(sys)");

    // Check our sysetm
    sys.checkInit();

    // Stop old solverthread
    stop();
    if(wait) 
      doSolve(this,&sys);
    else {
      stop();
      solverThread.reset(new std::thread(doSolve,this,&sys));
    }
  }
  void Solver::solve(TaSystem& sys,const SolverConfiguration& sc1){
    TRACE(20,"Solver::solve(sys,sc)");
    // Update solver configuration
    SolverConfiguration::operator=(sc1);
    solve(sys);
  }
  ErrorVals doIter(TaSystem* sys1,SolverConfiguration* sc) {
    TRACE(15,"Solver::doIter(sc)");
    using arma::norm;
    if(!sys1)
      throw MyError("No TaSystem given!");
    // In case no solverconfig is given
    std::unique_ptr<SolverConfiguration> sc_local;
    if(!sc){
      sc_local.reset(new SolverConfiguration());
      sc=sc_local.get();
    }
      
    // Wrap to reference
    TaSystem& sys=(*sys1);
    // old error
    vd error=sys.Error();
    // old error norm
    d oldfuner=norm(error);
    if(!(error.size()>0)){
      throw MyError("Error illegal residual vector obtained. Exiting.");
    }
    vd oldx=sys.getRes();

    us Ndofs=error.size();
    sp_mat jac=sys.jac(sc->dampfac);

    assert(jac.n_cols==error.size());
    assert(jac.n_rows==error.size());

    TRACE(15,"Solving linear system...");
    vd fulldx; 
    // This can throw. Is catched a layer higher
    fulldx=-1.0*arma::spsolve(jac,error,"superlu");

    d newfuner;
    d reler;
    vd newx(Ndofs);

    reler=norm(fulldx);

    do{
      newx=oldx+sc->dampfac*fulldx;
      sys.setRes(newx);
      newfuner=norm(sys.Error());
      if((newfuner>oldfuner || !(newfuner>0)) && sc->dampfac>sc->mindampfac){
        if(sc->dampfac*0.5<sc->mindampfac){
          d temporary=sc->mindampfac;
          sc->dampfac=temporary;
        }
        else
          sc->dampfac=sc->dampfac*0.5; // *= is not working (is atomic variable)
        cout << "Decreasing dampfac, new dampfac = " << sc->dampfac << "\n";
      }
    } while((newfuner>oldfuner || !(newfuner>0)) && sc->dampfac>sc->mindampfac);
    if(newfuner<oldfuner && sc->dampfac<sc->maxdampfac){
      cout << "Increasing dampfac, new dampfac = " << sc->dampfac << " . Max dampfac: "<< sc->maxdampfac << "\n";
      sc->dampfac=sc->dampfac*2;
    }

    else{
      newx=oldx+sc->dampfac*fulldx;
      sys.setRes(newx);
      newfuner=norm(sys.Error());
    }
    
    cout << "Current dampfac: " << sc->dampfac << "\n";
    TRACE(10,"Iteration done...");
    
    return {newfuner,reler};		// Return function error

  } // Solver::DoIter()

  Solver::~Solver()  {
    TRACE(25,"Solver::~Solver()");
    stop();
  }
} // namespace tasystem
