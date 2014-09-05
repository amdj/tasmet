#include "solver.h"
#include "solverinstance.h"
#include "solverconfiguration.h"
#include "vtypes.h"
#include <Eigen/Sparse>
#include "arma_eigen.h"

#include <boost/thread.hpp>


namespace tasystem{
  
  


  // A solver always contains a valid system.
  Solver::Solver(const TaSystem& sys):tasystem(sys.copy()) {
    TRACE(15,"Solver(TaSystem&)");
  }
  Solver::Solver(const Solver& o): Solver(o.sys()){}
  Solver& Solver::operator=(const Solver& other){
    tasystem.reset(other.sys().copy());
    return *this;
  }
 

  typedef tuple<d,d> dtuple;
  void Solver::solve(us maxiter,d funtol,d reltol,d dampfac){
    TRACE(20,"Solver started.");
    SolverConfiguration sc(maxiter,funtol,reltol,dampfac);
    SolverInstance si(*this,sc);
    boost::thread solverThread(si);
    // This results in segfault
    // solverThread.detach();
    // si.Start();

    // while(si.isRunning()){
    //   cout << "Solver runnning...\n";
    //   sleep(1);

    // }
    // si.Join();
    solverThread.join();
    cout << "Solver done.\n";
  }
  tuple<d,d> Solver::doIter(d dampfac){
    // Do an iteration
    assert(dampfac>0 && dampfac<=1.0);
    TRACE(15,"Solver::DoIter()");
    TRACE(10,"Computing error...");
    evd error=sys().error();
    assert(error.size()>0);
    TRACE(10,"Getting old result vector...");
    evd oldx=sys().getRes();
    // cout << "oldx:" << oldx <<"\n";
    d funer=error.norm();
    us Ndofs=error.size();
    TRACE(15,"Computing Jacobian...");
    esdmat jac=sys().jac();
    jac.makeCompressed();

    assert(jac.cols()==error.size());
    assert(jac.rows()==error.size());
    TRACE(15,"Solving linear system...");
    tuple<d,d> progres;
    try{
      evd dx=-1.0*dampfac*math_common::solvesys_eigen(jac,error);
      d reler=dx.norm();
      TRACE(10,"Setting new solution vector...");
      evd newx=oldx+dx;
      using math_common::armaView;
      vd Newx=armaView(newx);
      sys().setRes(Newx);
      progres=std::make_tuple(funer,reler);
    }
    catch(int e){
      throw e;
    }
    // cout << "dx:"<<dx<<"\n";
    TRACE(10,"Iteration done...");
    return progres;		// Return function error
  } // Solver::DoIter()

} // namespace tasystem
