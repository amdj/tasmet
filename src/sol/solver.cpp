#include "solver.h"
#include "tasystem.h"
#include "solverinstance.h"
#include "vtypes.h"
#include <Eigen/Sparse>


namespace tasystem{
  
  typedef std::tuple<d,d> dtuple;  
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
  std::tuple<d,d> Solver::doIter(d dampfac){
    TRACE(15,"Solver::doIter("<<dampfac<<")");
    if(dampfac>0 && dampfac<=1.0) // if dampfac is legal value, use that
      sc.dampfac=dampfac;


    evd error=sys().error();
    if(!(error.size()>0)){
      WARN("Error illegal residual vector obtained. Exiting.");
      return std::make_tuple(-1,-1);
    }
    evd oldx=sys().getRes();
    d oldfuner=error.norm();

    us Ndofs=error.size();

    esdmat jac;

    jac=sys().jac(sc.dampfac);
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
    TRACE(15,"SFSG");
    d newfuner;
    d reler;
    evd newx(Ndofs);
    vd Newx=math_common::armaView(newx);
    TRACE(15,"SFSG");
    reler=fulldx.norm();
    do{
      newx=oldx+sc.dampfac*fulldx;
      sys().setRes(Newx);
      newfuner=sys().error().norm();
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
      sys().setRes(Newx);
      newfuner=sys().error().norm();
    }
    
    cout << "Current dampfac: " << sc.dampfac << "\n";
    TRACE(10,"Iteration done...");
    return std::make_tuple(newfuner,reler);		// Return function error

  } // Solver::DoIter()

  Solver::~Solver()  {
    TRACE(15,"Solver::~Solver()");
    TRACE(25,"Waiting for solver to stop..");
    stop();
    delete tasystem;
  }
} // namespace tasystem
