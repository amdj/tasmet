#include "enginesystem.h"

namespace tasystem{



  
  EngineSystem::EngineSystem(const Globalconf& gc,TimingConstraint tc):
    TaSystem(gc),
    tc(tc)
  {    TRACE(15,"EngineSystem::EngineSystem(gc,tc)");}
  EngineSystem::EngineSystem(const EngineSystem& sys):
    TaSystem(sys),
    tc(sys.tc)
  {    TRACE(15,"EngineSystem::EngineSystem(EngineSystem))");}
  EngineSystem& EngineSystem::operator=(const EngineSystem& sys){
    TRACE(15,"EngineSystem::operator=()");
    TaSystem::operator=(sys);
    tc=sys.tc;
    
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
    mass*=gc.rho0;
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
    us segdofs;
    us startdof=0;
    us Nsegs=getNSegs();
    TRACE(-1,"Nsegs:"<< Nsegs);
    for(us i=0;i<Nsegs;i++){
      segdofs=segs[i]->getNDofs();
      dmtotdx.subvec(startdof,startdof+segdofs-1)=segs[i]->dmtotdx();
      startdof=startdof+segdofs;
    }
    return dmtotdx;
  }
  esdmat EngineSystem::jac(){
    TRACE(15,"Enginesystem::jac()");
    us Ndofs;
    esdmat subjac=TaSystem::jac();
    vtriplet jactriplets=math_common::getTripletsBlock(subjac,1,0,subjac.rows()-1,subjac.cols());

    if(gc.Nf>0)
      Ndofs=getNDofs()+1;
    else
      Ndofs=getNDofs();
    esdmat jac(Ndofs,Ndofs);

    
    math_common::shiftTriplets(jactriplets,1,0);
    // Reserve some extra space in this tripletlist and add the extra components
    vd dmtotdx=this->dmtotdx();
    // cout << "dmtotdx:" <<dmtotdx;
    us dmtotdxsize=dmtotdx.size();
    us newspace=jactriplets.size();

    newspace+=getNVertex();

    if(gc.Nf>0){
      vd domg=this->domg();
      // TRACE(50,"domg:"<<domg);

      us domgsize=domg.size();
      newspace+=domgsize;	// A little bit too much, but anyway..
      jactriplets.reserve(newspace);
      for(us i=0;i<domgsize;i++){
	if(domg(i)!=0)
	  jactriplets.push_back(triplet(i,Ndofs-1,domg(i)));
      }	// for domg

      // Now we need to add the timingconstraint as well.
      jactriplets.push_back(triplet(Ndofs-1,tc.dofnr(*this),1));

      
    } // gc.Nf>0
    jactriplets.reserve(newspace);
    for(us i=0;i<dmtotdxsize;i++)
      if(dmtotdx(i)!=0)
	jactriplets.push_back(triplet(0,i,dmtotdx(i)));
    
    jac.setFromTriplets(jactriplets.begin(),jactriplets.end());
    // cout << "Last column of jac: " << jac.innerVector(jac.cols()-1) <<"\n";
    return jac;
  }
  evd EngineSystem::error(){
    TRACE(15,"EngineSystem::error()");
    us Ndofs;
    if(gc.Nf>0)
      Ndofs=getNDofs()+1;
    else
      Ndofs=getNDofs();    
    evd error(Ndofs);		// Add one for the timing constraint
    if(gc.Nf>0)    {
      error.head(Ndofs-1)=TaSystem::error();
      error(Ndofs-1)=tc.error(*this);
    } else{
      error=TaSystem::error();
    }
    // Strip first equation (for now, assuming it is a continuity equation!!)
    error(0)=getCurrentMass()-gc.getMass();
    return error;
  }
  evd EngineSystem::getRes(){
    TRACE(15,"EngineSystem::getRes()");
    checkInit();
    if(gc.Nf>0){
      us Ndofs=getNDofs()+1;
      evd res(Ndofs);		// Add one for the timing constraint

      res.head(Ndofs-1)=TaSystem::getRes();
      res(Ndofs-1)=gc.getomg();
      return res;
    } else{
      return TaSystem::getRes();
    }
  }
  void EngineSystem::setRes(const vd& res){
    TRACE(15,"EngineSystem::setRes()");
    checkInit();
    // Sanity check
    us ndofs=getNDofs();
    us ressize=res.size();
    TRACE(50,"ndofs: "<< ndofs);
    TRACE(50,"res size: "<< ressize);
    assert(ressize==ndofs || ressize==ndofs+1);

    
    if(ressize==ndofs)
      TaSystem::setRes(res);
    else{
      assert(res(ndofs)>0);
      TaSystem::setRes(res.subvec(0,ndofs-1));
      gc.setomg(res(ndofs));
      TRACE(18,"New freq:"<< res(ndofs)/2/number_pi);
    }
  }
  
  vd EngineSystem::domg(){
    TRACE(15,"EngineSystem::domg()");
    assert(segs.size());
    us ndofs=getNDofs();
    vd domg(ndofs,fillwith::zeros);
    us segdofs;
    us startdof=0;
    us Nsegs=getNSegs();
    for(us i=0;i<Nsegs;i++){
      segdofs=segs[i]->getNDofs();
      domg.subvec(startdof,startdof+segdofs-1)=segs[i]->domg();
      startdof=startdof+segdofs;
    }
    // Eigen::SparseVector()
    return domg;
  }

} // namespace tasystem

