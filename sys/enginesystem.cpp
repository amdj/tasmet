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
  esdmat EngineSystem::jac(){
    TRACE(15,"Enginesystem::jac()");
    us Ndofs=getNDofs()+1;
    esdmat jac=TaSystem::jac();
    jac.resize(Ndofs,Ndofs);
    
    return jac;
  }
  evd EngineSystem::error(){
    TRACE(15,"EngineSystem::error()");
    if(gc.Nf>0){
      us Ndofs=getNDofs()+1;
      evd error(Ndofs);		// Add one for the timing constraint
      error.head(Ndofs-1)=TaSystem::error();
      error(Ndofs-1)=tc.error(*this);
      return error;
    } else{
      return TaSystem::error();
    }
  }
  evd EngineSystem::getRes(){
    TRACE(15,"EngineSystem::getRes()");
    if(gc.Nf>0){
      us Ndofs=getNDofs()+1;
      evd res(Ndofs);		// Add one for the timing constraint
      res.head(Ndofs-1)=TaSystem::getRes();
      res(Ndofs-1)=gc.omg;
      return res;
    } else{
      return TaSystem::getRes();
    }
  }


} // namespace tasystem

