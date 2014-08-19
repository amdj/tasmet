#include "solver.h"
#include "vtypes.h"
#include <Eigen/Sparse>
#include "arma_eigen.h"

namespace tasystem{
  using math_common::armaView;
  using Eigen::ComputationInfo;
  
  evd solvesys_eigen(const esdmat& K,const evd& f)  {
    TRACE(15,"Eigen solver used for large system");
    // Solve a linear system using Eigen sparse
    // Form: K*x=f
    assert(f.size()>0);

    // Eigen::SimplicialCholesky<esdmat,Eigen::COLAMDOrdering<int> > solver(jac2); // Werkt niet...
    // Eigen::SimplicialCholesky<esdmat,3 > solver(jac2); // Werkt niet...    
    // Eigen::SimplicialLDLT<esdmat> solver(jac2);
    // Eigen::SparseQR<esdmat,Eigen::COLAMDOrdering<int> > solver(jac2);      
    // TRACE(10,"Converting data to Eigen...");
    // esdmat jac2=sys.Jac();
    TRACE(10,"Creating solver...");
    // esdmat eig_K=math_common::ArmaToEigen(sdmat(K)); // Eigen
    // matrix

    // cout << "Matrix k:"<< K << "\n";
    
    Eigen::SparseLU<esdmat> solver(K);
    switch(solver.info()){
    case ComputationInfo::InvalidInput:
      cout << "Solver initialization failed: invalid input" << "\n";
      exit(1);
    case ComputationInfo::NumericalIssue:
      cout << "Solver initialization failed: numerical issue" << "\n";
      exit(1);
    case ComputationInfo::NoConvergence:
      cout << "Solver initialization failed: no convergence" << "\n";
      exit(1);

      
    }

    // evd f=math_common::ArmaToEigen(f);

    TRACE(10,"Solving linear system...");
    evd x=solver.solve(f);
    return x;    
  }
  
  // A solver always contains a valid system.
  Solver::Solver(const taSystem& sys):tasystem(sys.copy()) {
    TRACE(15,"Solver(taSystem&)");
  }
  Solver::Solver(const Solver& o): Solver(o.sys()){}
  Solver& Solver::operator=(const Solver& other){
    tasystem.reset(other.sys().copy());
    return *this;
  }
 

  typedef tuple<d,d> dtuple;
  void Solver::solve(us maxiter,d funtol,d reltol){
    TRACE(20,"Solver started.");

    evd&& error=sys().error();
    d funer=error.norm();
    d reler=0;
    us nloop=0;
    if(maxiter==0)
      maxiter=SOLVER_MAXITER;
    TRACE(15,"maxiter:"<< maxiter);
    TRACE(15,"funtol:"<< funtol);
    TRACE(15,"reltol:"<< reltol);
    while((funer>funtol || reler>reltol) && nloop<maxiter)
      {
	dtuple ers=doIter();
	funer=std::get<0>(ers);
	reler=std::get<1>(ers);
	cout << "Iteration: "<<nloop<<" , function error: "<<funer<<" , relative error:" << reler<< ".\n";
	nloop++;
      }
    if(nloop==maxiter)
      WARN("Solver reached maximum number of iterations! Results might not be reliable!");
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
    d funer=error.norm();
    us Ndofs=error.size();
    TRACE(15,"Computing Jacobian...");
    esdmat jac=sys().jac();
    jac.makeCompressed();

    assert(jac.cols()==error.size());
    assert(jac.rows()==error.size());
    TRACE(15,"Solving linear system...");
    evd dx=-1.0*dampfac*solvesys_eigen(jac,error);
    d reler=dx.norm();
    
    TRACE(10,"Setting new solution vector...");
    evd newx=oldx+dx;
    vd Newx=armaView(newx);

    cout << "dx:"<<dx<<"\n";
    sys().setRes(Newx);
    TRACE(10,"Iteration done...");
    return std::make_tuple(funer,reler);		// Return function error
  } // Solver::DoIter()

} // namespace tasystem
