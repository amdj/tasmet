#include <cassert>
#include "tasystem.h"
#include "triplets.h"
#include "jacobian.h"
#include "seg.h"
#include "connector.h"
#include "tube.h"
#include "utils.h"

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

  Seg* TaSystem::getSeg(us i) const { return (*this)[i];}
  
  Seg* TaSystem::operator[](us i) const {
    us nSegs=this->nSegs();
    if(i<nSegs)
      return segs[i];
    else
      return nullptr;
  }
  void TaSystem::setDofEqNrs(){
    TRACE(14,"TaSystem::setDofnrs()");
    us firstdof=0;
    us firsteq=0;

    for(auto seg=segs.begin();seg!=segs.end();seg++) {
      // Set the dofnrs
      VARTRACE(10,firstdof);
      VARTRACE(10,firsteq);
      (*seg)->setDofNrs(firstdof);
      (*seg)->setEqNrs(firsteq);
      firstdof+=(*seg)->getNDofs();
      firsteq+=(*seg)->getNEqs();
    }
    for(auto con=connectors.begin();con!=connectors.end();con++) {
      VARTRACE(10,firsteq)
      (*con)->setEqNrs(firsteq);
      firsteq+=(*con)->getNEqs();
    }
  }
  void TaSystem::checkInit(){		// Often called simple method: inline
    if(!hasInit){
      init();
    }
  }
  void TaSystem::init(){
    TRACE(24,"TaSystem::init()");
    hasInit=false;
    us Nsegs=nSegs();
    // TRACE(25,"Address gc:" <<&gc);
    us i=0;

    // Quite some assumptions are done where the order of this
    // initialization depends on. So first the segs, then the connectors.
    for(auto seg=segs.begin();seg!=segs.end();seg++) {
      // (*seg)->init(*this);
    }
    for(auto con=connectors.begin();con!=connectors.end();con++) {
      // (*con)->init(*this);
    }

    // Set all dofs and equation numbers
    setDofEqNrs();

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
  void TaSystem::setNf(us Nf){
    checkInit();
    TRACE(30,"TaSystem::setNf()");
    gc_.setNf(Nf);
    for(auto seg=segs.begin();seg!=segs.end();seg++)      {
      (*seg)->updateNf();
    }
    for(auto con=connectors.begin();con!=connectors.end();con++)      {
      (*con)->updateNf();
    }    TRACE(25,"New Ns:"<< gc_.Ns() << " . Adres gc: " << &gc_);

    setDofEqNrs();
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
  TripletList TaSystem::jacTriplets() {
    TRACE(14,"TaSystem::jacTriplets()");
    Jacobian j(getNDofs());
    us Nsegs=nSegs();
    for(const Seg* seg :segs)
      seg->jac(j);
    for(auto con: connectors)
      con->jac(j);

    return j;
  }
  sp_mat TaSystem::jac(d dummy){
    TRACE(14,"TaSystem::Jac()");
    checkInit();
    us Ndofs=getNDofs();
    TripletList triplets=jacTriplets();
    // Remove all invalid elements
    triplets.setValid();
    return triplets;            // Implicitly converts to sp_mat
  }
  vd TaSystem::Error() {
    TRACE(14,"TaSystem::Error()");
    checkInit();

    us Ndofs=getNDofs();
    vd Error(Ndofs); // Globally, neqs=ndofs
    Error.zeros();
    us Nsegs=nSegs();
    us Ncon=nConnectors();
    us segeqs,coneqs;
    us starteq=0;
    TRACE(-1,"Nsegs:"<< Nsegs);
    for(us i=0;i<Nsegs;i++){
      segeqs=segs[i]->getNEqs();
      Error.subvec(starteq,starteq+segeqs-1)=segs[i]->error();
      starteq+=segeqs;
    }
    for(us i=0;i<Ncon;i++){
      TRACE(15,"Connector " << i);
      coneqs=connectors[i]->getNEqs();
      Error.subvec(starteq,starteq+coneqs-1)=connectors[i]->error();
      starteq+=coneqs;
    }
    return Error;
  }
  vd TaSystem::getRes(){
    TRACE(14,"TaSystem::getRes()");
    checkInit();
    us Ndofs=getNDofs();
    TRACE(14,"TaSystem::GetRes(), Ndofs:"<< Ndofs);
    const us& Ns=gc_.Ns();
      
    vd Res(Ndofs);

    us segdofs;
    us startdof=0;
    us Nsegs=nSegs();
    for(us i=0;i<Nsegs;i++){
      segdofs=segs[i]->getNDofs();
      Res.subvec(startdof,startdof+segdofs-1)=segs[i]->getRes();
      startdof+=segdofs;
      TRACE(4,"Seg:"<<i<<", Ndofs: "<<segdofs);
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
  void TaSystem::setRes(const vd& Res){
    checkInit();
    TRACE(14,"TaSystem::SetRes(vd res)");
    us Ndofs=getNDofs();

    if(Res.size()==Ndofs){
      us segdofs;
      us startdof=0;
      us Nsegs=nSegs();
      for(us i=0;i<Nsegs;i++){
        segdofs=segs[i]->getNDofs();
        segs[i]->setRes(Res.subvec(startdof,startdof+segdofs-1));
        startdof=startdof+segdofs;
      }
    }
    else {
      WARN("Amount of DOFS in result vector does not match system size!\nNothing done.");
    }
  } // TaSystem::SetRes()
  d TaSystem::getCurrentMass() {
    TRACE(10,"TaSystem::getCurrentMass()");
    assert(hasInit);
    d mass=0;
    us nsegs=segs.size();
    for(us i=0;i<nsegs;i++){
      mass+=segs[i]->getCurrentMass();
    } // for loop
    return mass;
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

