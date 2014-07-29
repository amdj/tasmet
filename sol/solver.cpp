#include "solver.h"
#include <vtypes.h>

namespace tasystem{

  vd solvesys_eigen(const dmat& K,const vd& f)  {
    // Solve a linear system using Eigen sparse
    // Form: K*x=f
    assert(f.size()>0);
    using math_common::esdmat;
    using math_common::evd;
    vd x(f.size())    ;
    try{
      // Eigen::SimplicialCholesky<esdmat,Eigen::COLAMDOrdering<int> > solver(jac2); // Werkt niet...
      // Eigen::SimplicialCholesky<esdmat,3 > solver(jac2); // Werkt niet...    
      // Eigen::SimplicialLDLT<esdmat> solver(jac2);
      // Eigen::SparseQR<esdmat,Eigen::COLAMDOrdering<int> > solver(jac2);      
      // TRACE(10,"Converting data to Eigen...");
      // esdmat jac2=sys->Jac();
      TRACE(10,"Creating solver...");
      esdmat eig_K=math_common::ArmaToEigen(sdmat(K)); // Eigen matrix
      Eigen::SparseLU<esdmat> solver(eig_K);
      evd eig_f=math_common::ArmaToEigen(f);

      TRACE(10,"Solving linear system...");
      evd edx=solver.solve(eig_f);
      x=math_common::EigenToArma(edx);
    }
    catch(...)
      {
	// Todo something useful
      }

    return x;    
  }
  

  Solver::Solver(const TAsystem& sys1) {
    TRACE(15,"Solver(TAsystem&)");
    sys=new TAsystem(sys1);
  }
  Solver::Solver(const Solver& other){
    if(other.sys!=NULL)
      sys=new TAsystem(*other.sys);
    else
      sys=NULL;
  }
  Solver& Solver::operator=(const Solver& other){
    if(sys!=NULL)
      delete sys;
    if(other.sys!=NULL)
      sys=new TAsystem(*other.sys);
    else
      sys=NULL;
    return *this;
  }
  
  Solver::~Solver(){
      TRACE(-5,"~Solver()");
      if(sys!=NULL){
	delete sys;
	sys=NULL;
      }
  }
  typedef tuple<d,d> dtuple;
  void Solver::solve(d funtol,d reltol){
    TRACE(20,"Solver started.");
    vd error=sys->Error();
    d funer=arma::norm(error,2);
    d reler=0;
    us nloop=0;
    while((funer>funtol || reler>reltol) && nloop<SOLVER_MAXITER)
      {
	dtuple ers=DoIter();
	funer=std::get<0>(ers);
	reler=std::get<1>(ers);
	cout << "Iteration: "<<nloop<<" , function error: "<<funer<<" , relative error:" << reler<< ".\n";
	nloop++;
      }
    if(nloop==SOLVER_MAXITER)
      WARN("Solver reached maximum number of iterations! Results might not be reliable!");
    cout << "Solver done.\n";
  }
  tuple<d,d> Solver::DoIter(d dampfac){
    // Do an iteration
    assert(dampfac>0 && dampfac<=1.0);
    TRACE(15,"Solver::DoIter()");
    TRACE(10,"Computing error...");
    vd err=sys->Error();
    assert(err.size()>0);
    TRACE(10,"Updating result vector...");
    vd oldx=sys->GetRes();
    
    d funer=arma::norm(err,2);
    us Ndofs=err.size();

    TRACE(15,"Computing Jacobian...");
    dmat jac(sys->Jac());
    assert(jac.n_cols==err.size());
    assert(jac.n_rows==err.size());
    vd dx;
    if(Ndofs>500){
      TRACE(10,"Using eigen sparse solver to solve system.");
      dx=-1.0*dampfac*solvesys_eigen(jac,err);
    }
    else{
      TRACE(10,"Using dense Arma solver to solve system.");
      dx=-1.0*dampfac*arma::solve(jac,err);
    }
    d reler=arma::norm(dx,2);
    TRACE(10,"Solving linear system done.");      
    sys->SetRes(oldx+dx);
    TRACE(10,"Iteration done...");
    return std::make_tuple(funer,reler);		// Return function error
  } // Solver::DoIter()

} // namespace tasystem
