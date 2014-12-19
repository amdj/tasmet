#include "tasystem.h"
#include "triplets.h"
#include "jacobian.h"
#include "seg.h"
#include "connector.h"

namespace tasystem{
  using segment::Seg;
  using segment::Connector;
  using math_common::armaView;

  inline us max(us s,us t){  return s? s>=t : t;}

  TaSystem::TaSystem(const Globalconf& gc):gc(gc){
    TRACE(14,"TaSystem::TaSystem(gc)");
    this->gc.setDriven(true);
  }
  TaSystem::TaSystem(const TaSystem& o)
  {
    TRACE(14,"TaSystem::TaSystem(TaSystem&)");
    copyTaSystem(o);
  }
  void TaSystem::setGc(const Globalconf& gc){
    TRACE(14,"TaSystem::setGc()");
    this->gc=gc;
    hasInit=false;
  }
  void TaSystem::copyTaSystem(const TaSystem& o){
    TRACE(14,"TaSystem::copyTaSystem()");
    cleanup();
    gc=o.gc;
    assert(nSegs()==0);
    for(us i=0;i<o.nSegs();i++) {
      TRACE(14,"Copying segment "<<i << "...");
      assert(o.getSeg(i)!=NULL);
      addSeg(*o.getSeg(i));
      assert(nSegs()==i+1);
    }
    for(us i=0;i<o.nConnectors();i++) {
      TRACE(14,"Copying connector "<<i << "...");
      assert(o.connectors[i]);
      addConnector(*o.connectors[i]);
      assert(nConnectors()==i+1);
    }

    // segConnections=o.segConnections;
    hasInit=false;

  }
  TaSystem& TaSystem::operator=(const TaSystem& other){
    TRACE(14,"TaSystem::operator=()");
    cleanup();
    copyTaSystem(other);
    return *this;
  }
  void TaSystem::cleanup(){
    for (us i=0; i < segs.size(); ++i) {
      delete segs[i];
    }
    for (us i=0; i < connectors.size(); ++i) {
      delete connectors[i];
    }
    segs.clear();
    connectors.clear();

    hasInit=false;
  }
  void TaSystem::addSeg(const std::vector<Seg*>& segs){
    TRACE(14,"TaSystem::addSeg()");
    for(auto seg=segs.begin();seg!=segs.end();seg++){
      if(*seg!=NULL) 
        addSeg(**seg);
    }
  }
  void TaSystem::addSeg(const Seg& seg){
    TRACE(14,"TaSystem::addseg()");
    hasInit=false;
    segs.emplace_back(seg.copy());
    segs[nSegs()-1]->setNumber(nSegs()-1);
  }
  void TaSystem::addConnector(const Connector& con){
    TRACE(14,"TaSystem::addConnector()");
    hasInit=false;
    connectors.emplace_back(con.copy());
    connectors[nSegs()-1]->setNumber(nSegs()-1);
  }

  Seg* TaSystem::getSeg(us i) const { return (*this)[i];}
  
  Seg* TaSystem::operator[](us i) const {
    us nSegs=this->nSegs();
    if(i<nSegs)
      return segs[i];
    else
      return NULL;
  }
  void TaSystem::setDofEqNrs(){
    TRACE(14,"TaSystem::setDofnrs()");
    us firstdof=0;
    us firsteq=0;

    for(auto seg=segs.begin();seg!=segs.end();seg++) {
      // Set the dofnrs
      (*seg)->setDofNrs(firstdof);
      (*seg)->setEqNrs(firsteq);
      firstdof+=(*seg)->getNDofs();
      firsteq+=(*seg)->getNEqs();
    }
    for(auto con=connectors.begin();con!=connectors.end();con++) {
      (*con)->setEqNrs(firsteq);
      firsteq+=(*con)->getNEqs();
    }

  }
  void TaSystem::init(){
    TRACE(14,"TaSystem::init()");
    us Nsegs=nSegs();
    // TRACE(25,"Address gc:" <<&gc);
    us i=0;
    WARN("No seg connections");

    for(auto seg=segs.begin();seg!=segs.end();seg++) {
      (*seg)->init(*this);
    }
    for(auto con=connectors.begin();con!=connectors.end();con++) {
      (*con)->init(*this);
    }

    // Do some post-sanity checks
    setDofEqNrs();
    us Ndofs=getNDofs();
    gc.setSys(this);
    TRACE(10,"Segment initialization done. Total NDofs:"<< Ndofs);
    if(Ndofs>MAXNDOFS)      {
      WARN("Way too many DOFS required: Ndofs=" <<Ndofs << ". Initialization failed\n");
      return;
    }

    if(getNDofs()!=getNEqs()){
      WARN("Ndofs on TaSystem level not equal to number of equations! Initialization failed.");
      WARN("Ndofs="<< getNDofs());
      WARN("Neqs ="<< getNEqs());
      return;
      // exit(1);
    }
    // Last, but not least: initialize a pointer to this tasystem in
    // globalconf
    hasInit=true;
  }
  void TaSystem::setNf(us Nf){
    TRACE(30,"TaSystem::setNf()");
    gc.setNf(Nf);
    for(auto seg=segs.begin();seg!=segs.end();seg++)      {
      (*seg)->updateNf();
    }
    TRACE(25,"New Ns:"<< gc.Ns() << " . Adres gc: " << &gc);

    setDofEqNrs();
  }

  
  us TaSystem::getNDofs() const  {
    TRACE(14,"TaSystem::getNDofs()");
    us Ndofs=0;
    for(us i=0;i<nSegs();i++)
      Ndofs+=segs.at(i)->getNDofs();
    return Ndofs;
  }
  us TaSystem::getNEqs() const  {
    TRACE(14,"TaSystem::getNDofs()");
    us Neqs=0;
    for(us i=0;i<nSegs();i++)
      Neqs+=segs.at(i)->getNEqs();
    for(us i=0;i<nConnectors();i++)
      Neqs+=connectors.at(i)->getNEqs();

    return Neqs;
  }
  void TaSystem::show(us detailnr){
    TRACE(10,"TaSystem::show()");
    checkInit();
    TRACE(14,"checkInit() done");
    cout << "########################## Showing TaSystem...\n";
    cout << "Showing Global configuration...\n";
    gc.show();
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
    const us& Ns=gc.Ns();
    us Ndofs=getNDofs();
    TRACE(-1,"Ndofs:"<<Ndofs);    
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
      coneqs=connectors[i]->getNEqs();
      Error.subvec(starteq,starteq+coneqs-1)=connectors[i]->error();
    }
    return error;
  }
  evd TaSystem::getRes(){
    TRACE(14,"TaSystem::getRes()");
    checkInit();
    us Ndofs=getNDofs();
    TRACE(14,"TaSystem::GetRes(), Ndofs:"<< Ndofs);
    const us& Ns=gc.Ns();
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
  void TaSystem::setRes(const TaSystem& other){
    TRACE(25,"TaSystem::setRes(TaSystem)");
    WARN("This only works for Tube segments so far");
    checkInit();
    us nsegs=nSegs();
    assert(other.nSegs()==nsegs);
    WARN("Not yet available, testing should be done in ")
    for(us i=0;i<nsegs;i++) {
      getSeg(i)->setRes(*other.getSeg(i));
    }
    gc=other.gc;
    
  }
  void TaSystem::resetHarmonics(){
    assert(!segs.empty());
    for(auto seg=segs.begin();seg!=segs.end();seg++)
      (*seg)->resetHarmonics();
  }
  void TaSystem::setRes(const evd& res){
    TRACE(15,"EngineSystem::setRes()");
    vd res2=math_common::armaView(res);
    setRes(res2);
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
    else
      {
        WARN("Amount of DOFS in result vector does not match system size!\nNothing done.");
      }
  } // TaSystem::SetRes()
  d TaSystem::getCurrentMass() {
    TRACE(10,"TaSystem::getCurrentMass()");
    checkInit();
    d mass=0;
    us nsegs=segs.size();
    for(us i=0;i<nsegs;i++){
      mass+=segs[i]->getCurrentMass();
    } // for loop
    return mass;
  }
  void TaSystem::showJac(bool force){
    TRACE(15,"TaSystem::showJac()");
    checkInit();
    esdmat jac=this->jac();
    if(force || jac.cols()<50){
      Eigen::IOFormat CleanFmt(2,0," ",";\n","","","[","]");
      Eigen::MatrixXd jacd(jac);
      cout << "Jacobian: \n" << jacd.format(CleanFmt) << "\n";
      cout << "Determinant of jacobian:" << jacd.determinant() << "\n";
    }
    else if(jac.cols()<1000){
      Eigen::MatrixXd jacd(jac);
      cout << "Jacobian too large to show, but determinant of jacobian is:" << jacd.determinant() << "\n";      
    }
    else 
      cout << "Jacobian size is too large to show.\n";
    cout << "Number of columns in Jacobian: " << jac.cols() << "\n";
    cout << "Number of rows    in Jacobian: " << jac.rows() << "\n";    
  }
  TaSystem::~TaSystem() {
    TRACE(-5,"~TaSystem()");
    cleanup();
  }

} // namespace tasystem
