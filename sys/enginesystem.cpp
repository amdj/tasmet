#include "enginesystem.h"

namespace tasystem{



  
  EngineSystem::EngineSystem(const Globalconf& gc):
    TaSystem(gc)
  {    TRACE(15,"EngineSystem::EngineSystem(gc,tc)");}
  EngineSystem::EngineSystem(const EngineSystem& sys):
    TaSystem(sys),
    tc(sys.tc),av(sys.av)
  {    TRACE(15,"EngineSystem::EngineSystem(EngineSystem))");}
  EngineSystem& EngineSystem::operator=(const EngineSystem& sys){
    TRACE(15,"EngineSystem::operator=()");
    TaSystem::operator=(sys);
    tc=sys.tc;
    av=sys.av;
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
  vtriplet EngineSystem::Ljac(){
    TRACE(15,"Enginesystem::Ljac()");

    esdmat subjac=TaSystem::jac();
    vtriplet jactriplets=math_common::getTripletsBlock(subjac,1,0,subjac.rows()-1,subjac.cols());

    us Ndofs=getNDofs();
    if(gc.Nf>0)
      Ndofs++;

    math_common::shiftTriplets(jactriplets,1,0);
    // Reserve some extra space in this tripletlist and add the extra components
    vd dmtotdx=this->dmtotdx();
    // cout << "dmtotdx:" <<dmtotdx;
    us dmtotdxsize=dmtotdx.size();
    us extraspace=getNVertex();

    if(gc.Nf>0){
      vd domg=this->domg();
      // TRACE(50,"domg:"<<domg);

      us domgsize=domg.size();
      extraspace+=domgsize;	// A little bit too much, but anyway..
      math_common::reserveExtraDofs(jactriplets,extraspace);
      for(us i=0;i<domgsize;i++){
	if(domg(i)!=0)
	  jactriplets.push_back(triplet(i,Ndofs-1,domg(i)));
      }	// for domg

      // Now we need to add the timingconstraint as well.
      // TRACE(50,"Dofnr:"<<tc.dofnr(*this));
      jactriplets.push_back(triplet(Ndofs-1,tc.dofnr(*this),1));
      
    } // gc.Nf>0
    else			// Only dofs for dmtotdx
      math_common::reserveExtraDofs(jactriplets,extraspace);

    for(us i=0;i<dmtotdxsize;i++)
      if(dmtotdx(i)!=0)
	jactriplets.push_back(triplet(0,i,dmtotdx(i)));

    return jactriplets;
  }
  esdmat EngineSystem::jac(){
    TRACE(15,"EngineSystem::jac()");
    
    vtriplet Mjac=this->Ljac(); // Its called Mjac, but here it is still Ljac
    d aval=av.value(*this);
    assert(aval!=0);		// Otherwise, something is wrong.
    math_common::multiplyTriplets(Mjac,1/aval); // Now its Mjac


    evd M=this->error();
    us valdof=av.dofnr(*this);
    math_common::reserveExtraDofs(Mjac,M.size());
    // If we add extra triplets to this vector, summation is done
    // according to the eigen documentation.
    // TRACE(50,"SFSG")
    for(us i=0;i<M.size();i++)
      if(M(i)!=0)
    	Mjac.push_back(triplet(i,valdof,-M(i)/aval));
    us Ndofs=getNDofs();	// This number is without extra omega dof
    if(gc.Nf>0)
      Ndofs+=1;
    esdmat jac(Ndofs,Ndofs);
    jac.setFromTriplets(Mjac.begin(),Mjac.end());
    return jac;
  }


  
  evd EngineSystem::error(){
    TRACE(15,"EngineSystem::error()");
    if(gc.Nf>0)    {
      us Ndofs=getNDofs()+1;
      evd error(Ndofs);		// Add one for the timing constraint
      error.head(Ndofs-1)=TaSystem::error();
      error(Ndofs-1)=tc.value(*this);
      // Strip first equation (for now, assuming it is a continuity
      d aval=av.value(*this);
      TRACE(50,"aval:"<<aval);
      error(0)=getCurrentMass()-gc.getMass();
      error*=(1/aval);		// Divide L by amplitude value to
				// avoid zero amplitude as solution
      return error;
    } else{
      return TaSystem::error();
    }

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
  void EngineSystem::setTimingConstraint(us segnr,us vertexnr,us varnr,us freqnr){
    TRACE(18,"EngineSystem::setTimingConstraint()");
    tc.set(segnr,vertexnr,varnr,freqnr);
  }
  void EngineSystem::setAmplitudeDof(us segnr,us vertexnr,us varnr,us freqnr)   {
    TRACE(18,"EngineSystem::setAmplitudeDof()");
    av.set(segnr,vertexnr,varnr,freqnr);
  }

} // namespace tasystem

