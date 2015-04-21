// #define TRACERPLUS (10)
#include "enginesystem.h"
#include "triplets.h"


// To set divide by amplitude on
// #define DIVAMPL
#define TIMINGCONSTRAINT
#define DOMGFAC (1e0)

#define MASSEQ (0)
// #define MASSEQ (0)

namespace tasystem{
  using arma::sp_mat;

  EngineSystem::EngineSystem(const Globalconf& gc):
    TaSystem(gc)
  {
    TRACE(15,"EngineSystem::EngineSystem(gc,tc)");
    setDriven(false);
  }
  EngineSystem::EngineSystem(const EngineSystem& sys):
    TaSystem(sys),
    tc(sys.tc),av(sys.av)
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
    // Only set initial mass when it is not already set by user by
    // using setInitialMass()
    if(getMass()==0)
      setMass(getInitialMassFromGc());
    hasInit=true;
  }
  d EngineSystem::getInitialMassFromGc() const {
    TRACE(15,"EngineSystem::getInitialMassFromGc()");
    d mass=0;
    us nsegs=segs.size();
    for(us i=0;i<nsegs;i++){
      mass+=segs[i]->getCurrentMass();
    } // for loop
    TRACE(20,"Volume of device: "<< mass<<" [m^3].");
    mass=mass*gc_.rho0();
    // TRACE(20,"Initial mass as computed from rho0="<<gc.rho0<<" [kg/m^3],\n and volume of device: "<< mass << " [kg].");
    return mass;
  }
  vd EngineSystem::dmtotdx() const{
    TRACE(15,"EngineSystem::dmtotdx()");
    // Should become a row vector, but anyway.
    vd dmtotdx(getNDofs(),fillwith::zeros);
    us Nsegs=nSegs();
    for(us i=0;i<Nsegs;i++){
      segs[i]->dmtotdx(dmtotdx);
      // WARN("dmtotdx: Not done function");
    }
    return dmtotdx;
  }
  vd EngineSystem::Error(){
    TRACE(20,"EngineSystem::Error()");

    checkInit();
    TRACE(15,"EngineSystem::error()");
    if(gc_.Nf()>0){
      d aval=av.value(*this);
      cout << "Current amplitude value: " << aval << "\n";
      cout << "Current frequency      : " << gc_.getfreq() << "\n";
    }

    #ifdef DIVAMPL
    return errorM();
    #else
    return errorL();
    #endif
  }
  sp_mat EngineSystem::jac(d dampfac){
    TRACE(15,"EngineSystem::jac("<< dampfac<<")");
    #ifdef DIVAMPL
    TripletList jactr=this->Mjac(dampfac);
    #else
    TripletList jactr=this->Ljac(dampfac);
    #endif

    us Ndofs=getNDofs();	// This number is without extra omega dof
    #ifdef TIMINGCONSTRAINT
    TRACE(15,"Timincontraint on!");
    if(gc_.Nf()>0)
      Ndofs++;
    #endif
    return jactr;
  }
  // INCLUDING DIVIDE BY AMPLITUDE


  TripletList EngineSystem::Ljac(d dampfac){
    TRACE(15,"Enginesystem::Ljac("<<dampfac<<")");
    us ndofs=getNDofs();
    TripletList jactriplets=jacTriplets();

    us Ndofs=getNDofs();
    #ifdef TIMINGCONSTRAINT
    if(gc_.Nf()>0)
      Ndofs++;
    #endif
    jactriplets.zeroOutRow(MASSEQ);  // Replace this equation with global
                                // mass conservation
    // Reserve some extra space in this tripletlist and add the extra components
    vd dmtotdx=this->dmtotdx();
    // cout << "dmtotdx:" <<dmtotdx;
    us dmtotdxsize=dmtotdx.size();
    us extraspace=0;
    #ifdef TIMINGCONSTRAINT
    if(gc_.Nf()>0){

      vd domg=this->domg();
      if(dampfac>0){
        TRACE(15,"Dampfac set:"<< dampfac << "\nNorm of domg:"<< arma::norm(domg));
        domg*=DOMGFAC/dampfac;
      }
      // TRACE(50,"domg:"<<domg);
      us domgsize=domg.size();
      extraspace+=domgsize;	// A little bit too much, but anyway..
      jactriplets.reserveExtraDofs(extraspace);
      for(us i=0;i<domgsize;i++){
        if(domg(i)!=0)
          jactriplets.push_back(Triplet(i,Ndofs-1,domg(i)));
      }	// for domg

      // Now we need to add the timingconstraint as well.
      TRACE(14,"Dofnr timingconstraint:"<<tc.dofnr(*this));
      jactriplets.push_back(Triplet(Ndofs-1,tc.dofnr(*this),1));
      
    } // gc.Nf()>0
    else			// Only dofs for dmtotdx
      #endif
      jactriplets.reserveExtraDofs(extraspace);

    for(us k=0;k<dmtotdxsize;k++)
      if(dmtotdx(k)!=0){
        // TRACE(20,"k: " << k);
        // TRACE(20,"dmtotdx:"<< dmtotdx(k));
        jactriplets.push_back(Triplet(MASSEQ,k,dmtotdx(k)));
      }
    return jactriplets;
  }
  TripletList EngineSystem::Mjac(d dampfac){
    TRACE(15,"EngineSystem::Mjac("<<dampfac<<")");
    
    TripletList Mjac=this->Ljac(dampfac); // Its called Mjac, but here it is still Ljac
    if(gc_.Nf()>0){
      d aval=av.value(*this);
      assert(aval!=0);		// Otherwise, something is wrong.
      Mjac.multiplyTriplets(1/aval); // Now its nearly Mjac

      vd M=this->errorM();
      us valdof=av.dofnr(*this);
      Mjac.reserveExtraDofs(M.size());
      // If we add extra triplets to this vector, summation is done
      // according to the eigen documentation.

      for(us i=0;i<M.size();i++)
        if(M(i)!=0)
          Mjac.push_back(Triplet(i,valdof,-M(i)/aval));
      // Now it is officially Mjac
    }
    return Mjac;
  }
  vd EngineSystem::errorM(){
    TRACE(15,"EngineSystem::errorM()");
    if(gc_.Nf()>0)    {
      d aval=av.value(*this);
      assert(aval!=0);
      // Divide L by amplitude value to
      // avoid zero amplitude as solution
      return errorL()*(1.0/aval);		// Add one for the timing constraint
    } else{
      return errorL();      
    }
  }
  // WITHOUT DIVIDE BY AMPLITUDE
  vd EngineSystem::errorL(){
    TRACE(15,"EngineSystem::errorL()");
    
    us Ndofs=getNDofs();
    #ifdef TIMINGCONSTRAINT
    if(gc_.Nf()>0)
      Ndofs++;
    #endif
    vd error(Ndofs);		// Add one for the timing constraint  
    #ifdef TIMINGCONSTRAINT
    if(gc_.Nf()>0)    {
      error.subvec(0,Ndofs-2)=TaSystem::Error();
      error(Ndofs-1)=tc.value(*this);
    }
    else {
      #endif
      error=TaSystem::Error();
    #ifdef TIMINGCONSTRAINT      
    }      
    #endif
    // Strip first equation (for now, assuming it is a continuity
    d curmass=getCurrentMass();
    TRACE(20,"Current mass in system: "<< curmass << " [kg] \n");
    error(MASSEQ)=curmass-getMass();
    return error;
  }
  vd EngineSystem::getRes(){
    TRACE(15,"EngineSystem::getRes()");
    checkInit();
    #ifdef TIMINGCONSTRAINT
    if(gc_.Nf()>0){
      us Ndofs=getNDofs()+1;
      vd res(Ndofs);		// Add one for the timing constraint

      res.subvec(0,Ndofs-2)=TaSystem::getRes();
      res(Ndofs-1)=gc_.getomg();
      return res;
    } else{
      #endif
      return TaSystem::getRes();
    #ifdef TIMINGCONSTRAINT
    }
    #endif
  }
  void EngineSystem::setRes(const vd& res){
    TRACE(15,"EngineSystem::setRes()");
    checkInit();
    // Sanity check
    us ndofs=getNDofs();
    us ressize=res.size();
    TRACE(15,"ndofs: "<< ndofs);
    TRACE(15,"res size: "<< ressize);
    assert(ressize==ndofs || ressize==ndofs+1);

    
    if(ressize==ndofs)
      TaSystem::setRes(res);
    else{
      TaSystem::setRes(res.subvec(0,ndofs-1));
      gc_.setomg(res(ndofs));
      TRACE(18,"New freq:"<< res(ndofs)/2/number_pi);
    }
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
  void EngineSystem::setTimingConstraint(us segnr,us cellnr,us Varnr,us freqnr){
    TRACE(18,"EngineSystem::setTimingConstraint()");
    tc.set(segnr,cellnr,Varnr,freqnr);
  }
  void EngineSystem::setAmplitudeDof(us segnr,us cellnr,us Varnr,us freqnr)   {
    TRACE(18,"EngineSystem::setAmplitudeDof()");
    av.set(segnr,cellnr,Varnr,freqnr);
  }

} // namespace tasystem

