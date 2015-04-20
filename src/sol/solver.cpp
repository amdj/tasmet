#include "solver.h"
#include "tasystem.h"
#include "solverinstance.h"
#include "vtypes.h"



namespace tasystem{
  
  typedef std::tuple<d,d> dtuple;  
  using arma::sp_mat;
  

  // A solver always contains a valid system.
  Solver::Solver(const TaSystem& sys):tasystem(sys.copy()) {
    TRACE(15,"Solver(TaSystem&)");
  }
  Solver::Solver(const Solver& o): Solver(*o.tasystem){}
  void Solver::stop() {
    TRACE(15,"Solver::stop()");
    if(solverThread){
      cout << "Waiting for solver to finish...\n";
      sc.maxiter=0;
      solverThread->join();
      solverThread.reset();
    }
  }


  void Solver::solve(us maxiter,d funtol,d reltol,d mindampfac,d maxdampfac,bool wait){
    TRACE(20,"Solver started.");
    sys().checkInit();
    sc=SolverConfiguration(maxiter,funtol,reltol,mindampfac,maxdampfac);

    // Stop old solverthread
    stop();

    solverThread.reset(new boost::thread(SolverInstance(*this)));

    if(wait){
      cout << "Waiting for solver...\n";
      solverThread->join();
      solverThread.reset();
    }

  }
  ErrorVals Solver::doIter(d dampfac){
    TRACE(15,"Solver::doIter("<<dampfac<<")");
    using arma::norm;
    if(dampfac>0 && dampfac<=1.0) // if dampfac is legal value, use that
      sc.dampfac=dampfac;


    vd error=sys().Error();
    if(!(error.size()>0)){
      WARN("Error illegal residual vector obtained. Exiting.");
      return {-1,-1};
    }
    vd oldx=sys().getRes();
    d oldfuner=norm(error);

    us Ndofs=error.size();

    sp_mat jac=sys().jac(sc.dampfac);

    // assert(jac.cols()==error.size());
    // assert(jac.rows()==error.size());

    TRACE(15,"Solving linear system...");
    vd fulldx; 
    try{
      
      fulldx=-1.0*arma::spsolve(jac,error,"superlu");
    }
    catch(int e){
      throw e;
    }
    TRACE(15,"SFSG");
    d newfuner;
    d reler;
    vd newx(Ndofs);

    reler=norm(fulldx);
    do{
      newx=oldx+sc.dampfac*fulldx;
      sys().setRes(newx);
      newfuner=norm(sys().Error());
      if((newfuner>oldfuner || !(newfuner>0)) && sc.dampfac>sc.mindampfac){
        sc.dampfac=sc.dampfac*0.5;
        cout << "Decreasing dampfac, new dampfac = " << sc.dampfac << "\n";
      }
    } while((newfuner>oldfuner || !(newfuner>0)) && sc.dampfac>sc.mindampfac);
    if(newfuner<oldfuner && sc.dampfac<sc.maxdampfac){
      cout << "Increasing dampfac, new dampfac = " << sc.dampfac << " . Max dampfac: "<< sc.maxdampfac << "\n";
      sc.dampfac=sc.dampfac*2;
    }

    else{
      newx=oldx+sc.dampfac*fulldx;
      sys().setRes(newx);
      newfuner=norm(sys().Error());
    }
    
    cout << "Current dampfac: " << sc.dampfac << "\n";
    TRACE(10,"Iteration done...");
    return {newfuner,reler};		// Return function error

  } // Solver::DoIter()

  Solver::~Solver()  {
    TRACE(15,"Solver::~Solver()");
    TRACE(25,"Waiting for solver to stop..");
    stop();
    delete tasystem;
  }
} // namespace tasystem
