// #define TRACERPLUS (10)
#include "exception.h"
#include "triplets.h"
#include "enginesystem.h"
#include "seg.h"
#include <cassert>

// To set divide by amplitude on
// #define DIVAMPL
#define TIMINGCONSTRAINT
#define DOMGFAC (1e0)

namespace tasystem{
  using arma::sp_mat;
  using segment::Seg;

  EngineSystem::EngineSystem(const Globalconf& gc):
    TaSystem(gc)
  {
    TRACE(15,"EngineSystem::EngineSystem(gc,tc)");
    // Since an EngineSystem cannot be driven:
    setDriven(false);
    if(gc.Nf()<1)
      throw MyError("Too low frequency number given. Should be at "
                    "least 1 for an EngineSystem");
  }
  EngineSystem::EngineSystem(const EngineSystem& sys):
    TaSystem(sys)
  {
    TRACE(15,"EngineSystem::EngineSystem(EngineSystem))");
    setDriven(false);
  }
  EngineSystem::EngineSystem(const TaSystem& sys):TaSystem(sys){}
  void EngineSystem::show(us detailnr) {
    checkInit();
    cout << "########################## Showing EngineSystem...\n";
    cout << "Target system mass: " << getMass() << "\n";
    TaSystem::show(detailnr);
  }
  void  EngineSystem::init(){
    TRACE(15,"EngineSystem::init()");
    TaSystem::init();
    hasInit=false;
    // Getting the PhaseConstraint Dof
    for(const Seg* seg: segs) {
      int thissegPhaseDof=seg->providePhaseDof();
      if(thissegPhaseDof>=0) {
        if(PhaseConstraintDof<0) {
          // set the PhaseConstraintDof to this seg its 
          PhaseConstraintDof=thissegPhaseDof;
          // Set the pointer to the right segment
          phaseDofSeg=seg;
        }
        else {
          throw MyError("Too many segments with a Phase Constraint"
                        " Dof! This can be maximally one segment.");
        }
      } // if thissegDhaseDof>0
    } // for

    // Check if we really found a PhaseConstraintDof
    if(PhaseConstraintDof<1){
      throw MyError("No phase constraints found. Initialization of "
                    "EngineSystem failed.");
    }

    // Crude check if the PhaseConstraintdof we found is valid
    if((us) PhaseConstraintDof>getNDofs()) {
        throw MyError("Invalid phase constraint Dof given."
                      " Is the frequency OK?");
    }

    // If everything is allright...
    hasInit=true;
  }
  vd EngineSystem::Error(){
    TRACE(20,"EngineSystem::Error()");
    checkInit();
    assert(phaseDofSeg);
    // Add one Dof for the timing constraint  
    us Ndofs=getNDofs()+1;

    vd error(Ndofs);
    error.subvec(0,Ndofs-1)=TaSystem::Error();
    error(Ndofs-1)=phaseDofSeg->PhaseDofValue();

    return error;
  }

  TripletList EngineSystem::jacTriplets(d dampfac){
    TRACE(15,"Enginesystem::Ljac("<<dampfac<<")");
    us ndofs=getNDofs()+1;
    // Get Jacobian from TaSystem
    TripletList jactriplets=TaSystem::jacTriplets(dampfac);
    // Increase number of Dofs with 1
    jactriplets.setNdofs(ndofs);
    us extraspace=0;
    vd domg=this->domg();
    us domgsize=domg.size();
    // if(dampfac>0){
    //   TRACE(15,"Dampfac set:"<< dampfac << "\nNorm of domg:"<< arma::norm(domg));
    //   domg*=DOMGFAC/dampfac;
    // }
    // TRACE(50,"domg:"<<domg);
    // extraspace+=domgsize;	// A little bit too much, but anyway..
    // jactriplets.reserveExtraDofs(extraspace);
    for(us i=0;i<domgsize;i++){
      if(domg(i)!=0)
        jactriplets.push_back(Triplet(i,ndofs-1,domg(i)));
    }	// for domg
    // Now we need to add the timingconstraint as well.
    TRACE(14,"Dofnr timingconstraint:"<<PhaseConstraintDof);
    jactriplets.push_back(Triplet(ndofs-1,PhaseConstraintDof,1));

    // jactriplets.reserveExtraDofs(extraspace);
    return jactriplets;
  }
  // TripletList EngineSystem::Mjac(d dampfac){
  //   TRACE(15,"EngineSystem::Mjac("<<dampfac<<")");
    
  //   TripletList Mjac=this->Ljac(dampfac); // Its called Mjac, but here it is still Ljac
  //   if(gc_.Nf()>0){
  //     d aval=av.value(*this);
  //     assert(aval!=0);		// Otherwise, something is wrong.
  //     Mjac.multiplyTriplets(1/aval); // Now its nearly Mjac

  //     vd M=this->errorM();
  //     us valdof=av.dofnr(*this);
  //     Mjac.reserveExtraDofs(M.size());
  //     // If we add extra triplets to this vector, summation is done
  //     // according to the eigen documentation.

  //     for(us i=0;i<M.size();i++)
  //       if(M(i)!=0)
  //         Mjac.push_back(Triplet(i,valdof,-M(i)/aval));
  //     // Now it is officially Mjac
  //   }
  //   return Mjac;
  // }
  // vd EngineSystem::errorM(){
  //   TRACE(15,"EngineSystem::errorM()");
  //   if(gc_.Nf()>0)    {
  //     d aval=av.value(*this);
  //     assert(aval!=0);
  //     // Divide L by amplitude value to
  //     // avoid zero amplitude as solution
  //     return errorL()*(1.0/aval);		// Add one for the timing constraint
  //   } else{
  //     return errorL();      
  //   }
  // }
  // WITHOUT DIVIDE BY AMPLITUDE
  vd EngineSystem::getRes(){
    TRACE(15,"EngineSystem::getRes()");
    checkInit();
    us Ndofs=getNDofs()+1;
    vd res(Ndofs);		// Add one for the timing constraint
    res.subvec(0,Ndofs-2)=TaSystem::getRes();
    res(Ndofs-1)=gc_.getomg();
    return res;
  }
  void EngineSystem::setRes(const vd& res){
    TRACE(15,"EngineSystem::setRes()");
    checkInit();
    // Sanity check
    us ndofs=getNDofs()+1;
    us ressize=res.size();
    if(ressize!=ndofs)
      throw MyError("Result vector has inappropriate length!");

    TaSystem::setRes(res.subvec(0,ndofs-2));
    gc_.setomg(res(ndofs-1));
    TRACE(18,"New freq:"<< res(ndofs)/2/number_pi);
  }
  vd EngineSystem::domg(){
    TRACE(15,"EngineSystem::domg()");
    checkInit();
    assert(segs.size());
    us ndofs=getNDofs();        // On system level, this is neqs
    vd domg(ndofs,fillwith::zeros);
    us segdofs;
    us startdof=0;
    us Nsegs=nSegs();
    for(us i=0;i<Nsegs;i++){
      segs[i]->domg(domg);
    }
    // Eigen::SparseVector()
    return domg;
  }

} // namespace tasystem

