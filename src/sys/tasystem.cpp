#include "tasystem.h"
#include "triplets.h"
#include "jacobian.h"
#include "seg.h"
#include "connector.h"
#include "tube.h"
#include "piston.h"
#include "utils.h"
#include <cassert>

namespace tasystem{
  using segment::Seg;
  using segment::Connector;
  using arma::sp_mat;


  TaSystem::TaSystem(const Globalconf& gc):gc_(gc){
    TRACE(14,"TaSystem::TaSystem(gc)");
    this->setDriven(true);
  }
  TaSystem::TaSystem(const TaSystem& o):
    gc_(o.gc_)
  {
    TRACE(25,"TaSystem::TaSystem(TaSystem&) copy");
    for(auto it=o.segs.begin();it!=o.segs.end();it++){
      (*this)+=*(*it);
    }
    for(auto it=o.connectors.begin();it!=o.connectors.end();it++){
      (*this)+=*(*it);
    }
    hasInit=false;
  }
  void TaSystem::setGc(const Globalconf& gc){
    TRACE(14,"TaSystem::setGc()");

    if(gc_.Nf()!=gc.Nf()){
      gc_.setNf(gc.Nf());
      updateNf(gc.Nf());
    }
    gc_=gc;
    
    hasInit=false;
  }
  void TaSystem::cleanup(){
    TRACE(25,"TaSystem::cleanup()");
    utils::purge(segs);
    utils::purge(connectors);
    hasInit=false;
  }
  TaSystem& TaSystem::operator+=(const Seg& seg){
    TRACE(14,"TaSystem::operator+=(Seg)");
    hasInit=false;
    segs.emplace_back(seg.copy(*this));
    segs[nSegs()-1]->setNumber(nSegs()-1);
    return *this;
  }
  TaSystem& TaSystem::operator+=(const Connector& con){
    TRACE(24,"TaSystem::operator+=(Connector)");
    hasInit=false;
    connectors.emplace_back(con.copy(*this));
    connectors[nConnectors()-1]->setNumber(nConnectors()-1);
    return *this;
  }

  const Seg* TaSystem::getSeg(us i) const { return segs.at(i);}
  const Connector* TaSystem::getConnector(us i) const { return connectors.at(i);}

  void TaSystem::init(){
    TRACE(14,"TaSystem::init()");
    hasInit=false;
    arbitrateMassEq=-1;
    us firstdof=0;
    us firsteq=0;
    d mass=0;
    for(auto seg=segs.begin();seg!=segs.end();seg++) {
      // Set the dofnrs
      VARTRACE(10,firstdof);
      VARTRACE(10,firsteq);
      (*seg)->setDofNrs(firstdof);
      (*seg)->setEqNrs(firsteq);
      firstdof+=(*seg)->getNDofs();
      firsteq+=(*seg)->getNEqs();

      mass+=(*seg)->getMass();

    }
    // Update the mass in the system if it is undefined   (set by the
    // initial invalid state of -1)
    if(mass_<0)
      setMass(mass);

    for(auto con=connectors.begin();con!=connectors.end();con++) {
      VARTRACE(10,firsteq)
      (*con)->setEqNrs(firsteq);
      firsteq+=(*con)->getNEqs();
    }

    // Iterate over all segments and connectors to find one segment of
    // Dof that is willing to arbitrate the mass in the system. This
    // can, of course, be only be one segment or connector.
    using segment::SegConBase;
    // Build an fill a vector of all segments and connectors
    vector<const SegConBase*> segcons;
    segcons.insert(segcons.end(),segs.begin(),segs.end());
    segcons.insert(segcons.end(),connectors.begin(),connectors.end());
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
      WARN("Ndofs="<< getNDofs());
      WARN("Neqs ="<< getNEqs());
      throw MyError("Ndofs on TaSystem level not equal to number of equations! Initialization failed.");
    }
    hasInit=true;
  }
  void TaSystem::checkInit(){		// Often called simple method: inline
    if(!hasInit){
      init();
    }
  }
  void TaSystem::updateNf(us Nf){
    checkInit();
    TRACE(30,"TaSystem::setNf()");
    gc_.setNf(Nf);
    for(auto seg=segs.begin();seg!=segs.end();seg++)      {
      (*seg)->updateNf();
    }
    for(auto con=connectors.begin();con!=connectors.end();con++)      {
      (*con)->updateNf();
    }
    TRACE(25,"New Ns:"<< gc_.Ns() << " . Adres gc: " << &gc_);

    init();
  }
  
  us TaSystem::getNDofs() const  {
    TRACE(0,"TaSystem::getNDofs()");
    us Ndofs=0;
    for (auto seg = segs.begin(); seg != segs.end(); ++seg) {
      Ndofs+=(*seg)->getNDofs();
    }  
    return Ndofs;
  }
  us TaSystem::getNEqs() const  {
    TRACE(0,"TaSystem::getNDofs()");
    us Neqs=0;

    for (auto seg = segs.begin(); seg != segs.end(); ++seg) {
      Neqs+=(*seg)->getNEqs();
    }
    for (auto con = connectors.begin(); con != connectors.end(); ++con) {
      Neqs+=(*con)->getNEqs();
    }
    return Neqs;
  }
  void TaSystem::show(us detailnr){
    TRACE(25,"TaSystem::show()");
    checkInit();
    cout << "########################## Showing TaSystem...\n";
    cout << "Showing Global configuration...\n";
    gc_.show();
    if(detailnr>0){
      cout << "Now showing connectors in TaSystem...\n";
      for(us i=0;i<nConnectors();i++){
        TRACE(13,"Showing connector"<<i <<"..");
        connectors[i]->show(detailnr);
      }
      cout << "Now showing segments in TaSystem...\n";
      for(us i=0;i<nSegs();i++){
        TRACE(13,"Showing segment "<<i <<"..");
        segs[i]->show(detailnr);
      }
    } // detailnr>0
  }
  TripletList TaSystem::jacTriplets(d dummy) {
    TRACE(14,"TaSystem::jacTriplets()");
    us ndofs=getNDofs();
    Jacobian j(ndofs);

    for(const Seg* seg :segs)
      seg->jac(j);
    for(auto con: connectors)
      con->jac(j);
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
    for(us i=0;i<Nsegs;i++){
      segeqs=segs[i]->getNEqs();
      error.subvec(starteq,starteq+segeqs-1)=segs[i]->error();
      starteq+=segeqs;
    }
    for(us i=0;i<Ncon;i++){
      TRACE(15,"Connector " << i);
      coneqs=connectors[i]->getNEqs();
      error.subvec(starteq,starteq+coneqs-1)=connectors[i]->error();
      starteq+=coneqs;
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
    for(const Seg* seg: segs){
      segdofs=seg->getNDofs();
      Res.subvec(startdof,startdof+segdofs-1)=seg->getRes();
      startdof+=segdofs;
    }
    return Res;
  }
  // void TaSystem::setRes(const TaSystem& other){
  //   TRACE(25,"TaSystem::setRes(TaSystem)");
  //   WARN("This only works for Tube segments so far");
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
    for(auto seg=segs.begin();seg!=segs.end();seg++)
      (*seg)->resetHarmonics();
  }
  const tube::Tube& TaSystem::getTube(us i)  const {
    return dynamic_cast<const tube::Tube&>(*segs.at(i));
  }
  const mech::Piston& TaSystem::getPiston(us i)  const {
    return dynamic_cast<const mech::Piston&>(*segs.at(i));
  }
  void TaSystem::setRes(const vd& Res){
    checkInit();
    TRACE(14,"TaSystem::SetRes(vd res)");
    us Ndofs=getNDofs();

    if(Res.size()==Ndofs){
      us segdofs;
      us startdof=0;
      us Nsegs=nSegs();
      for(Seg* seg: segs){
        segdofs=seg->getNDofs();
        seg->setRes(Res.subvec(startdof,startdof+segdofs-1));
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
    for(const Seg* seg: segs){
      mass+=seg->getMass();
    } // for loop
    return mass;
  }
  vd TaSystem::dmtotdx() const{
    TRACE(15,"EngineSystem::dmtotdx()");
    // Should become a row vector, but anyway.
    vd dmtotdx(getNDofs(),fillwith::zeros);
    us Nsegs=nSegs();
    for(const Seg* seg:segs){
      seg->dmtotdx(dmtotdx);
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

