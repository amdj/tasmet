// #define TRACERPLUS (10)
#include "enginesystem.h"
#include "triplets.h"


// To set divide by amplitude on
#define DIVAMPL
#define TIMINGCONSTRAINT
#define DOMGFAC (1e-1)

#define MASSEQ (0)
// #define MASSEQ (0)

namespace tasystem{
  EngineSystem::EngineSystem(const Globalconf& gc):
    TaSystem(gc)
  {
    TRACE(15,"EngineSystem::EngineSystem(gc,tc)");
    this->gc.setDriven(false);
  }
  EngineSystem::EngineSystem(const EngineSystem& sys):
    TaSystem(sys),
    tc(sys.tc),av(sys.av)
  {
    TRACE(15,"EngineSystem::EngineSystem(EngineSystem))");
    gc.setDriven(false);
  }
  EngineSystem::EngineSystem(const TaSystem& sys):TaSystem(sys){}
  EngineSystem& EngineSystem::operator=(const EngineSystem& sys){
    TRACE(15,"EngineSystem::operator=()");
    TaSystem::operator=(sys);
    tc=sys.tc;
    av=sys.av;
    return *this;
  }
  void EngineSystem::show(bool showvertices) {
    checkInit();
    cout << "########################## Showing EngineSystem...\n";
    cout << "Target system mass: " << gc.getMass() << "\n";

    TaSystem::show(showvertices);

  }
  void  EngineSystem::init(){
    TRACE(15,"EngineSystem::init()");
    TaSystem::init();
    // Only set initial mass when it is not already set by user by
    // using setInitialMass()
    if(gc.getMass()==0)
      gc.setMass(getInitialMassFromGc());
  }
  d EngineSystem::getInitialMassFromGc() const {
    TRACE(15,"EngineSystem::getInitialMassFromGc()");
    d mass=0;
    us nsegs=segs.size();
    for(us i=0;i<nsegs;i++){
      mass+=segs[i]->geom.getFluidVolume();
    } // for loop
    TRACE(20,"Volume of device: "<< mass<<" [m^3].");
    mass*=gc.rho0;
    TRACE(20,"Initial mass as computed from rho0="<<gc.rho0<<" [kg/m^3],\n and volume of device: "<< mass << " [kg].");
    return mass;
  }
  us EngineSystem::getNVertex() const{
    us nvertex=0;
    for(auto seg=segs.begin();seg!=segs.end();++seg)
      nvertex+=seg->get()->getNVertex();
    return nvertex;
  }
  vd EngineSystem::dmtotdx() const{
    TRACE(15,"EngineSystem::dmtotdx()");
    // Should become a row vector, but anyway.
    vd dmtotdx(getNDofs(),fillwith::zeros);
    us Nsegs=getNSegs();
    for(us i=0;i<Nsegs;i++){
      segs[i]->dmtotdx(dmtotdx);
      // WARN("dmtotdx: Not done function");
    }
    return dmtotdx;
  }
  evd EngineSystem::error(){
    checkInit();
    TRACE(15,"EngineSystem::error()");
    if(gc.Nf()>0){
      d aval=av.value(*this);
      cout << "Current amplitude value: " << aval << "\n";
      cout << "Current frequency      : " << gc.getfreq() << "\n";
    }

    #ifdef DIVAMPL
    return errorM();
    #else
    return errorL();
    #endif
  }
  esdmat EngineSystem::jac(d dampfac){
    TRACE(15,"EngineSystem::jac("<< dampfac<<")");
    #ifdef DIVAMPL
    TripletList jactr=this->Mjac(dampfac);
    #else
    TripletList jactr=this->Ljac(dampfac);
    #endif

    us Ndofs=getNDofs();	// This number is without extra omega dof
    #ifdef TIMINGCONSTRAINT
    TRACE(15,"Timincontraint on!");
    if(gc.Nf()>0)
      Ndofs++;
    #endif
    esdmat jac(Ndofs,Ndofs);

    jactr.setValid();
    jac.setFromTriplets(jactr.begin(),jactr.end());
    return jac;
  }
  // INCLUDING DIVIDE BY AMPLITUDE


  TripletList EngineSystem::Ljac(d dampfac){
    TRACE(15,"Enginesystem::Ljac("<<dampfac<<")");

    TripletList jactriplets;
    jacTriplets(jactriplets);

    us Ndofs=getNDofs();
    #ifdef TIMINGCONSTRAINT
    if(gc.Nf()>0)
      Ndofs++;
    #endif
    jactriplets.zeroOutRow(MASSEQ);  // Replace this equation with global
                                // mass conservation
    // Reserve some extra space in this tripletlist and add the extra components
    vd dmtotdx=this->dmtotdx();
    // cout << "dmtotdx:" <<dmtotdx;
    us dmtotdxsize=dmtotdx.size();
    us extraspace=getNVertex();
    #ifdef TIMINGCONSTRAINT
    if(gc.Nf()>0){

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
    if(gc.Nf()>0){
      d aval=av.value(*this);
      assert(aval!=0);		// Otherwise, something is wrong.
      Mjac.multiplyTriplets(1/aval); // Now its nearly Mjac

      evd M=this->errorM();
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
  evd EngineSystem::errorM(){
    TRACE(15,"EngineSystem::errorM()");
    if(gc.Nf()>0)    {
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
  evd EngineSystem::errorL(){
    TRACE(15,"EngineSystem::errorL()");
    
    us Ndofs=getNDofs();
    #ifdef TIMINGCONSTRAINT
    if(gc.Nf()>0)
      Ndofs++;
    #endif
    evd error(Ndofs);		// Add one for the timing constraint  
    #ifdef TIMINGCONSTRAINT
    if(gc.Nf()>0)    {
      error.head(Ndofs-1)=TaSystem::error();
      error(Ndofs-1)=tc.value(*this);
    }
    else {
      #endif
      error=TaSystem::error();
    #ifdef TIMINGCONSTRAINT      
    }      
    #endif
    // Strip first equation (for now, assuming it is a continuity
    d curmass=getCurrentMass();
    TRACE(20,"Current mass in system: "<< curmass << " [kg] \n");
    error(MASSEQ)=curmass-gc.getMass();
    return error;
  }
  evd EngineSystem::getRes(){
    TRACE(15,"EngineSystem::getRes()");
    checkInit();
    #ifdef TIMINGCONSTRAINT
    if(gc.Nf()>0){
      us Ndofs=getNDofs()+1;
      evd res(Ndofs);		// Add one for the timing constraint

      res.head(Ndofs-1)=TaSystem::getRes();
      res(Ndofs-1)=gc.getomg();
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
      gc.setomg(res(ndofs));
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
    us Nsegs=getNSegs();
    for(us i=0;i<Nsegs;i++){
      segs[i]->domg(domg);
    }
    // Eigen::SparseVector()
    return domg;
  }
  void EngineSystem::setTimingConstraint(us segnr,us vertexnr,us varnr,us freqnr){
    TRACE(18,"EngineSystem::setTimingConstraint()");
    tc.set(segnr,vertexnr,varnr,freqnr);
  }
  void EngineSystem::setAmplitudeDof(us segnr,us vertexnr,us varnr,us freqnr)   {
    TRACE(18,"EngineSystem::setAmplitudeDof()");
    av.set(segnr,vertexnr,varnr,freqnr);
  }

} // namespace tasystem

