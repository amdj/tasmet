#include "solver.h"
#include "tasystem.h"
#include "solverinstance.h"
#include "vtypes.h"
#include <Eigen/Sparse>
#include "arma_eigen.h"


namespace tasystem{
  
  
  using Eigen::ComputationInfo;
  
  evd solvesys_eigen(const esdmat& K,const evd& f) {
    TRACE(15,"Eigen solver used for large system");
    // Solve a linear system using Eigen sparse
    // Form: K*x=f
    assert(f.size()>0);

    // Eigen::SimplicialCholesky<esdmat,Eigen::COLAMDOrdering<int> > solver(jac); // Werkt niet...
    // Eigen::SimplicialCholesky<esdmat,3 > solver(jac2); // Werkt niet...    
    // Eigen::SimplicialLDLT<esdmat> solver(jac2);
    // Eigen::SparseQR<esdmat,Eigen::COLAMDOrdering<int> > solver(jac2);      
    // matrix
    // cout << "f:\n"<<f;
    // cout << "Matrix k:"<< K << "\n";
    // Eigen::FullPivLU<edmat> dec(K);
    TRACE(19,"Initializing solver...");    
    // Eigen::SparseQR<esdmat,Eigen::COLAMDOrdering<int> > solver(K);
    Eigen::SparseLU<esdmat,Eigen::COLAMDOrdering<int> > solver(K);
    switch(solver.info()){
    case ComputationInfo::InvalidInput:
      cout << "Solver initialization failed: invalid input" << "\n";
      throw(0);
    case ComputationInfo::NumericalIssue:
      cout << "Solver initialization failed: numerical issue" << "\n";
      throw(0);
    case ComputationInfo::NoConvergence:
      cout << "Solver initialization failed: no convergence" << "\n";
      throw(0);
    case ComputationInfo::Success:
      break;
    } // switch
    TRACE(19,"Logarithm of absolute value of determinant of matrix: "<<solver.logAbsDeterminant());
    // Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver(K);
    // cout << "Logarithm of absolute value of determinant of matrix: "<<solver.logAbsDeterminant() <<"\n";
    
    // solver.setMaxIterations(500);
    
    TRACE(19,"Solving linear system...");    
    evd x=solver.solve(f);
    // std::cout << "#iterations: " << solver.iterations() << std::endl;
    // std::cout << "estimated error: " << solver.error() << std::endl;
 
    TRACE(19,"Solving linear system done.");    
    return x;    
  }


  // A solver always contains a valid system.
  Solver::Solver(const TaSystem& sys):tasystem(sys.copy()),dampfac(1.0) {
    TRACE(15,"Solver(TaSystem&)");
  }
  Solver::Solver(const Solver& o): Solver(*o.tasystem.get()){}
  Solver& Solver::operator=(const Solver& other){
    tasystem.reset(other.tasystem->copy());
    return *this;
  }
  void Solver::stop() {
    TRACE(15,"Solver::stop()");
    if(solverThread){
      cout << "Waiting for solver to finish...\n";
      sc->maxiter=0;
      solverThread->join();
      solverThread.reset();
    }
  }

  typedef tuple<d,d> dtuple;
  void Solver::solve(us maxiter,d funtol,d reltol,d mindampfac,bool wait){
    TRACE(20,"Solver started.");
    this->mindampfac=mindampfac;
    this->dampfac=1.0;
    if(maxiter==0)
      maxiter=SOLVER_MAXITER;

    if(solverThread){
      if(solverThread->joinable()){
        cout << "Killing old solver thread...\n";
        stop();
      }
    }
    
    if(!solverThread){
      sc.reset(new SolverConfiguration(maxiter,funtol,reltol,dampfac));
      // SolverInstance si(*this,*sc);
      solverThread.reset(new boost::thread(SolverInstance(*this,*sc)));
      // si();
      // This results in segfault
      // solverThread.detach();
      // si.Start();

      if(wait){
        cout << "Waiting for solver...\n";
        solverThread->join();
        solverThread.reset();
      }
      // si.Join();
      // solverThread.join();
    } else{
      WARN("Solver already running!");
    }

  }
  tuple<d,d> Solver::doIter(d dampfac1){
    TRACE(15,"Solver::doIter("<<dampfac1<<")");

    // Do an iteration
    if(dampfac1>0)
      this->dampfac=dampfac1;


    assert(dampfac>0 && dampfac<=1.0);

    evd error=sys().error();
    assert(error.size()>0);
    evd oldx=sys().getRes();
    d oldfuner=error.norm();

    us Ndofs=error.size();
    TRACE(15,"Computing Jacobian...");
    esdmat jac=sys().jac(dampfac);
    jac.makeCompressed();

    assert(jac.cols()==error.size());
    assert(jac.rows()==error.size());

    TRACE(15,"Solving linear system...");
    evd fulldx; 
    try{
      fulldx=-1.0*solvesys_eigen(jac,error);
    }
    catch(int e){
      throw e;
    }

    d newfuner;
    d reler;
    evd newx(Ndofs);
    vd Newx=math_common::armaView(newx);
    reler=fulldx.norm();
    do{
      TRACE(25,"Dampfac:"<< dampfac);
      newx=oldx+dampfac*fulldx;
      sys().setRes(Newx);
      newfuner=sys().error().norm();
      if((newfuner>oldfuner && dampfac>mindampfac)){
        dampfac=dampfac*0.9;
        cout << "Decreasing dampfac, new dampfac = " << dampfac << "\n";
      }
    } while((dampfac>mindampfac && newfuner>oldfuner));
    if(dampfac<=0.91 && newfuner<oldfuner)
      dampfac=dampfac/0.91;
    cout << "Current dampfac: " << dampfac << "\n";
    TRACE(10,"Iteration done...");
    return std::make_tuple(newfuner,reler);		// Return function error

  } // Solver::DoIter()

  Solver::~Solver()  {
    stop();
  }
} // namespace tasystem
