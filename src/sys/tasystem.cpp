#include "tasystem.h"
#include "triplets.h"
#include "jacobian.h"
#include "seg.h"
#include "connector.h"
#include "connectorvolume.h"
#include "duct.h"
#include "piston.h"
#include "utils.h"
#include <cassert>
#include <stdexcept>
#include "staticmsg.h"

namespace tasystem{
  using segment::Seg;
  using segment::Connector;
  using arma::sp_mat;
  common::StaticMsg<> ermsg;

  TaSystem::TaSystem(const Globalconf& gc):Globalconf(Globalconf::airSTP(0,100)) {
    TRACE(14,"TaSystem::TaSystem(gc)");
    Globalconf::operator=(gc);
  }
  TaSystem::TaSystem(const TaSystem& o):Globalconf(o)
  {
    TRACE(25,"TaSystem::TaSystem(TaSystem&) copy");
    hasInit=false;
    // First a constructor is called. After that, a function init() is
    // called. in TaSystem::init(). This way, virtual functions of the object can be called.

    for(auto it=o.segs.begin();it!=o.segs.end();it++){
      segs[it->first]=it->second->copy(*this);
    }
    for(auto it=o.connectors.begin();it!=o.connectors.end();it++){
      connectors[it->first]=it->second->copy(*this);
    }
  }
  void TaSystem::setGc(const Globalconf& gc){
    TRACE(14,"TaSystem::setGc()");

    if(Nf()!=gc.Nf()){
      updateNf(gc.Nf());
    }
    Globalconf::operator=(gc);
    
    init();
  }
  void TaSystem::cleanup(){
    TRACE(25,"TaSystem::cleanup()");
    utils::purge(segs);
    utils::purge(connectors);
    hasInit=false;
  }
  TaSystem& TaSystem::operator+=(const Seg& s){
    TRACE(24,"TaSystem::operator+=(Seg)");
    hasInit=false;
    segs[s.getID()]=s.copy(*this);
    segs[s.getID()]->setNumber(nSegs()-1);
    return *this;
  }
  TaSystem& TaSystem::operator+=(const Connector& c){
    TRACE(24,"TaSystem::operator+=(Connector)");
    hasInit=false;
    connectors[c.getID()]=c.copy(*this);
    connectors[c.getID()]->setNumber(nConnectors()-1);
    return *this;
  }

  const Seg* TaSystem::getSeg(const string& i) const { return segs.at(i);}
  const Connector* TaSystem::getConnector(const string& i) const { return connectors.at(i);}

  void TaSystem::init(){
    TRACE(14,"TaSystem::init()");
    hasInit=false;
    arbitrateMassEq=-1;
    us firstdof=0;
    us firsteq=0;
    d mass=0;
    for(auto seg: segs) {
      // Set the dofnrs
      seg.second->init();
      VARTRACE(10,firstdof);
      VARTRACE(10,firsteq);
      (seg.second)->setDofNrs(firstdof);
      (seg.second)->setEqNrs(firsteq);
      firstdof+=(seg.second)->getNDofs();
      firsteq+=(seg.second)->getNEqs();

      mass+=(seg.second)->getMass();

    }
    // Update the mass in the system if it is undefined   (set by the
    // initial invalid state of -1)
    if(mass_<0)
      setMass(mass);

    for(auto con: connectors) {
      con.second->init();
      con.second->setEqNrs(firsteq);
      firsteq+=con.second->getNEqs();
    }

    // Iterate over all segments and connectors to find one segment of
    // Dof that is willing to arbitrate the mass in the system. This
    // can, of course, be only be one segment or connector.
    using segment::SegConBase;
    // Build an fill a vector of all segments and connectors
    vector<const SegConBase*> segcons;
    for(auto seg:segs)
      segcons.push_back(seg.second);
    for(auto con:connectors)
      segcons.push_back(con.second);

    // Find a probable equation for mass arbitration. This should
    // indeed be done AFTER all segments and connectors have equation
    // and DOF numbers, as done in the previous for loops
    for(const SegConBase* segcon: segcons) {
      d segArbitrateMassEq=segcon->arbitrateMassEq();
      if(segArbitrateMassEq!=-1) {
        TRACE(15,"Found a segment of connector to arbitrate mass");
        if(arbitrateMassEq!=-1) {
          throw MyError("Found another segment that wants to arbitrate"
                        "the mass in the system. Only one segment or "
                        "connector can do that");
        }
        else {
          arbitrateMassEq=segArbitrateMassEq;
        }
      } // if segArbitratemassEq
        
    } // for segcon


    // Do some post-sanity checks
    us Ndofs=getNDofs();

    TRACE(10,"Segments initialization done. Total NDofs:"<< Ndofs);
    if(Ndofs>constants::maxndofs)      {
      throw MyError("Way too many DOFS required. Initialization failed.");
    }

    if(getNDofs()!=getNEqs()){
      throw MyError(ermsg("Ndofs on TaSystem level (%d) not equal to number "
                           "of equations (%d)! Initialization failed.",getNDofs(),getNEqs()));
    }
    hasInit=true;
  }
  void TaSystem::checkInit(){		// Often called simple method: inline
    if(!hasInit){
      init();
    }
  }
  void TaSystem::updateNf(us Nf){
    TRACE(30,"TaSystem::updateNf()");
    checkInit();
    setNf(Nf);
    for(auto seg: segs)      {
      seg.second->updateNf();
    }
    for(auto con: connectors)      {
      con.second->updateNf();
    }
    TRACE(25,"New Ns:"<< Ns());

    init();
  }
  
  us TaSystem::getNDofs() const  {
    TRACE(0,"TaSystem::getNDofs()");
    us Ndofs=0;
    for (auto seg : segs) {
      Ndofs+=seg.second->getNDofs();
    }  
    return Ndofs;
  }
  us TaSystem::getNEqs() const  {
    TRACE(0,"TaSystem::getNDofs()");
    us Neqs=0;

    for (auto seg :segs) {
      Neqs+=seg.second->getNEqs();
    }
    for (auto con :connectors) {
      Neqs+=con.second->getNEqs();
    }
    return Neqs;
  }
  void TaSystem::show(us detailnr){
    TRACE(25,"TaSystem::show()");
    checkInit();
    cout << "########################## Showing TaSystem...\n";
    cout << "Showing Global configuration...\n";
    Globalconf::show();
    if(detailnr>0){
      for(auto con:connectors){
	cout << "Showing connector with ID " << con.first << "\n";
        con.second->show(detailnr);
      }
      for(auto seg:segs){
	cout << "Showing segment with ID " << seg.first << "\n";
        seg.second->show(detailnr);
      }
    } // detailnr>0
  }
  TripletList TaSystem::jacTriplets(d dummy) {
    TRACE(14,"TaSystem::jacTriplets()");
    us ndofs=getNDofs();
    Jacobian j(ndofs);

    for(auto seg :segs)
      seg.second->jac(j);
    for(auto con: connectors)
      con.second->jac(j);
    // Convert to tripletlist
    TripletList jac=j;

    assert(arbitrateMassEq< (int) ndofs);
    // Exchange equation if we need to arbitrate mass
    if(arbitrateMassEq!=-1) {
      // Replace this equation with global mass conservation
      jac.zeroOutRow(arbitrateMassEq);
      vd dmtotdx=this->dmtotdx();
      us dmtotdxsize=dmtotdx.size();
      for(us k=0;k<dmtotdxsize;k++)
        if(dmtotdx(k)!=0){
          // TRACE(20,"k: " << k);
          // TRACE(20,"dmtotdx:"<< dmtotdx(k));
          jac.push_back(Triplet(arbitrateMassEq,k,dmtotdx(k)));
        }
    }
    return jac;
  }
  sp_mat TaSystem::jac(d dampfac){
    TRACE(14,"TaSystem::Jac()");
    checkInit();
    us Ndofs=getNDofs();
    TripletList triplets=jacTriplets(dampfac);
    // Remove all invalid elements
    triplets.setValid();
    return triplets;            // Implicitly converts to sp_mat
  }
  vd TaSystem::Error() {
    TRACE(14,"TaSystem::Error()");
    checkInit();

    us Ndofs=getNDofs();
    vd error(Ndofs); // Globally, neqs=ndofs
    error.zeros();
    us Nsegs=nSegs();
    us Ncon=nConnectors();
    us segeqs,coneqs;
    us starteq=0;
    TRACE(-1,"Nsegs:"<< Nsegs);
    us i=0;			// iterator
    for(auto seg:segs){
      TRACE(15,"Segment " << i);
      segeqs=seg.second->getNEqs();
      error.subvec(starteq,starteq+segeqs-1)=seg.second->error();
      starteq+=segeqs;
      i++;
    }
    i=0;
    for(auto con: connectors){
      TRACE(15,"Connector " << i);
      coneqs=con.second->getNEqs();
      error.subvec(starteq,starteq+coneqs-1)=con.second->error();
      starteq+=coneqs;
      i++;
    }

    assert(arbitrateMassEq< (int) Ndofs);
    // Exchange equation if we need to arbitrate mass
    if(arbitrateMassEq!=-1) {
      error(arbitrateMassEq)=getCurrentMass()-getMass();
    }

    return error;
  }
  vd TaSystem::getRes(){
    TRACE(14,"TaSystem::getRes()");
    checkInit();
    us Ndofs=getNDofs();
    vd Res(Ndofs);
    us startdof=0,segdofs;
    for(auto seg: segs){
      segdofs=seg.second->getNDofs();
      Res.subvec(startdof,startdof+segdofs-1)=seg.second->getRes();
      startdof+=segdofs;
    }
    return Res;
  }
  // void TaSystem::setRes(const TaSystem& other){
  //   TRACE(25,"TaSystem::setRes(TaSystem)");
  //   WARN("This only works for Duct segments so far");
  //   if(!checkInit())
  //     return;
  //   us nsegs=nSegs();
  //   assert(other.nSegs()==nsegs);
  //   WARN("Not yet available, testing should be done in ")
  //   for(us i=0;i<nsegs;i++) {
  //     getSeg(i)->setRes(*other.getSeg(i));
  //   }
  //   gc_=other.gc_;
  // }
  void TaSystem::resetHarmonics(){
    checkInit();
    assert(!segs.empty());
    for(auto seg: segs)
      seg.second->resetHarmonics();
  }
  const duct::Duct& TaSystem::getDuct(const string& id)  const {
    return dynamic_cast<const duct::Duct&>(*segs.at(id));
  }
  const duct::ConnectorVolume& TaSystem::getConnnectorVolume(const string& id)  const {
    return dynamic_cast<const duct::ConnectorVolume&>(*segs.at(id));
  }
  const mech::Piston& TaSystem::getPiston(const string& id)  const {
    return dynamic_cast<const mech::Piston&>(*segs.at(id));
  }
  void TaSystem::setRes(const vd& Res){
    checkInit();
    TRACE(14,"TaSystem::SetRes(vd res)");
    us Ndofs=getNDofs();

    if(Res.size()==Ndofs){
      us segdofs;
      us startdof=0;
      us Nsegs=nSegs();
      for(auto seg: segs){
        segdofs=seg.second->getNDofs();
        seg.second->setRes(Res.subvec(startdof,startdof+segdofs-1));
        startdof+=segdofs;
      } // for
      // Update the mass of the system from the result of the segments
      if(arbitrateMassEq==-1)
        setMass(getCurrentMass());
    }
    else {
      WARN("Amount of DOFS in result vector does not match system size!\nNothing done.");
    }
  } // TaSystem::SetRes()
  d TaSystem::getCurrentMass() {
    TRACE(10,"TaSystem::getCurrentMass()");
    checkInit();
    d mass=0;
    for(auto seg: segs){
      mass+=seg.second->getMass();
    } // for loop
    return mass;
  }
  vd TaSystem::dmtotdx() const{
    TRACE(15,"EngineSystem::dmtotdx()");
    // Should become a row vector, but anyway.
    vd dmtotdx(getNDofs(),fillwith::zeros);
    us Nsegs=nSegs();
    for(auto seg:segs){
      seg.second->dmtotdx(dmtotdx);
    }
    return dmtotdx;
  }

  dmat TaSystem::showJac(){
    TRACE(15,"TaSystem::showJac()");
    checkInit();
    return dmat(jac());
  }
  TaSystem::~TaSystem() {
    TRACE(25,"~TaSystem()");
    cleanup();
  }
} // namespace tasystem

