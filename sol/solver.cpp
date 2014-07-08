#include "solver.h"
#include <vtypes.h>

namespace tasystem{


  Solver::Solver(TAsystem& sys1): sys(sys1){

  }
  void Solver::DoIter(d dampfac){
    // Do an iteration
    using math_common::esdmat;
    using math_common::evd;
    TRACE(10,"Computing error...");
    vd err=sys->Error();
    TRACE(10,"Computing Jacobian...");
    sdmat jac(sys->Jac());
    // *jac=sdmat(sys->Jac());
    // dmat jac=Jac();    
    TRACE(10,"Converting data to Eigen...");
    // esdmat jac2=sys->Jac();
    esdmat jac2=math_common::ArmaToEigen(jac);
    // delete jac;
    evd Eerr=math_common::ArmaToEigen(err);

    // TRACE(10,"jac2 (eigen):"<<endl<<jac2);
    TRACE(10,"Computing old x...");
    vd oldx=sys->GetRes();
    try{
    // Eigen::SimplicialCholesky<esdmat,Eigen::COLAMDOrdering<int> > solver(jac2); // Werkt niet...
    // Eigen::SimplicialCholesky<esdmat,3 > solver(jac2); // Werkt niet...    
      TRACE(10,"Creating solver...");
      Eigen::SparseLU<esdmat> solver(jac2);
    // Eigen::SimplicialLDLT<esdmat> solver(jac2);
    
    // Eigen::SparseQR<esdmat,Eigen::COLAMDOrdering<int> > solver(jac2);      
      TRACE(10,"Solving linear system...");
      evd edx=solver.solve(Eerr);
      vd dx=-dampfac*math_common::EigenToArma(edx);

      // vd dx=-solve(dmat(jac),err);
      TRACE(10,"Solving linear done.");
      
      TRACE(10,"Updating result vector...");
      sys->SetRes(oldx+dx);
      TRACE(10,"Iteration done...");
    }
    catch(...)
      {
	
      }
  } // Solver::DoIter()

} // namespace tasystem
