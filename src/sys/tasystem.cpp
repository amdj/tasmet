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
  using math_common::armaView;

  inline us max(us s,us t){  return s? s>=t : t;}

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
  void TaSystem::jacTriplets(TripletList& trips){
    TRACE(14,"TaSystem::jacTriplets()");
    Jacobian jnew;
    us Nsegs=nSegs();
    for(auto seg=segs.begin();seg!=segs.end();seg++){
      (*seg)->jac(jnew);
    } // end for loop
    for(auto con=connectors.begin();con!=connectors.end();con++){
      (*con)->jac(jnew);
    } // end for loop

    // TRACE(25,"Jac\n"<<jac);
    trips=jnew.getTriplets();
    // trips.show();
    // cout << "Ndofs:" << getNDofs() << "\n";

  }
  esdmat TaSystem::jac(d dummy){
    TRACE(14,"TaSystem::Jac()");
    checkInit();

    // Something interesting has to be done here later on to connect
    // the different segments in the sense that blocks of Jacobian
    // matrix parts have to be moved to the right place etc. To be
    // continued...
    const us& Ns=gc_.Ns();
    us Ndofs=getNDofs();
    TRACE(15,"Ndofs:"<<Ndofs);    
    TripletList triplets;
    jacTriplets(triplets);

    triplets.setValid();
    esdmat eigsjac(Ndofs,Ndofs);
    eigsjac.setFromTriplets(triplets.begin(),triplets.end());
    return eigsjac;
  }
  evd TaSystem::error(){
    TRACE(14,"TaSystem::Error()");
    checkInit();

    us Ndofs=getNDofs();
    evd error(Ndofs);
    vd Error(error.data(),Ndofs,false,false); // Globally, neqs=ndofs
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
    // VARTRACE(15,error)
    return error;
  }
  evd TaSystem::getRes(){
    TRACE(14,"TaSystem::getRes()");
    checkInit();
    us Ndofs=getNDofs();
    TRACE(14,"TaSystem::GetRes(), Ndofs:"<< Ndofs);
    const us& Ns=gc_.Ns();
    evd res(Ndofs);
    
    vd Res(res.data(),Ndofs,false,false);

    us segdofs;
    us startdof=0;
    us Nsegs=nSegs();
    for(us i=0;i<Nsegs;i++){
      segdofs=segs[i]->getNDofs();
      Res.subvec(startdof,startdof+segdofs-1)=segs[i]->getRes();
      startdof+=segdofs;
      TRACE(4,"Seg:"<<i<<", Ndofs: "<<segdofs);
    }
    return res;
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
  void TaSystem::setRes(const evd& res){
    TRACE(15,"EngineSystem::setRes()");
    vd res2=math_common::armaView(res);
    setRes(res2);
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
    Eigen::MatrixXd jacd(jac());
    TRACE(15,"Eigen matrix built");
    return math_common::EigenToArma(jacd);
  }
  TaSystem::~TaSystem() {
    TRACE(25,"~TaSystem()");
    cleanup();
  }
  vd TaSystem::GetRes(){
    return math_common::EigenToArma(getRes());
  }
} // namespace tasystem

