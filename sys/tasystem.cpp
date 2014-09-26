#include "tasystem.h"
#include "triplets.h"

namespace tasystem{
  using segment::SegBase;
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
  void TaSystem::copyTaSystem(const TaSystem& o){
    TRACE(14,"TaSystem::copyTaSystem()");
    gc=Globalconf(o.gc);
    assert(getNSegs()==0);
    for(us i=0;i<o.getNSegs();i++)
      {
	TRACE(14,"Copying segment "<<i << "...");
	assert(o.getSeg(i)!=NULL);
	addSeg(*o.getSeg(i));
	assert(getNSegs()==i+1);
      }
    segConnections=o.segConnections;
    hasInit=false;

  }
  TaSystem& TaSystem::operator=(const TaSystem& other){
    TRACE(14,"TaSystem::operator=()");
    cleanup();
    copyTaSystem(other);
    return *this;
  }
  void TaSystem::cleanup(){
    segs.clear();
    segConnections.clear();
    hasInit=false;
  }
  void TaSystem::addSeg(const SegBase& seg){
    TRACE(14,"TaSystem::addseg()");
    hasInit=false;
    segs.emplace_back(seg.copy());
    segs[getNSegs()-1]->setNumber(getNSegs()-1);
  }
  SegBase* TaSystem::getSeg(us i) const { return (*this)[i];}
  
  SegBase* TaSystem::operator[](us i) const {
    us nSegs=getNSegs();
    if(i<nSegs)
      return segs[i].get();
    else
      return NULL;
  }

  void TaSystem::init(){
    TRACE(14,"TaSystem::init()");
    us Nsegs=getNSegs();
    us Ndofs=0;

    us i=0;
    for(auto v=segConnections.begin();v!=segConnections.end();++v){
      TRACE(90,"Connecting segment connection " << i << "...");
      coupleSegs(*v,*this);
      i++;
    }
    us firstdof=0;
    us firsteq=0;

    for(us i=0;i<Nsegs;i++)      {
      TRACE(9,"Initializing Segment "<<i<<"...");
      assert(segs.at(i));
      segs.at(i)->init(gc);

      // Set the dofnrs
      segs.at(i)->setDofNrs(firstdof);
      segs.at(i)->setEqNrs(firsteq);
      us thisndofs=segs.at(i)->getNDofs();
      us thisneqs=segs.at(i)->getNEqs();      
      Ndofs+=thisndofs;
      TRACE(12,"Ndofs for segment "<< i << ": "<<thisndofs);
      TRACE(12,"Neqs for segment "<< i << ": "<<thisneqs);      
      firstdof+=thisndofs;
      firsteq+=segs.at(i)->getNEqs();

    }

    // Do some post-sanity checks
    
    TRACE(10,"Segment initialization done. Total NDofs:"<< Ndofs);
    if(Ndofs>MAXNDOFS)      {
      WARN("Way too many DOFS required: Ndofs=" <<Ndofs << ". Exiting...\n");
      exit(1);
    }
    if(getNDofs()!=getNEqs()){
      WARN("Ndofs on TaSystem level not equal to number of equations!");
      WARN("Ndofs="<< getNDofs());
      WARN("Neqs ="<< getNEqs());
      exit(1);
    }
    
    // Last, but not least: initialize a pointer to this tasystem in
    // globalconf
    gc.setSys(this);
    
    hasInit=true;
  }
  

  
  us TaSystem::getNDofs() const  {
    TRACE(14,"TaSystem::getNDofs()");
    us Ndofs=0;
    for(us i=0;i<getNSegs();i++)
      Ndofs+=segs.at(i)->getNDofs();
    return Ndofs;
  }
  us TaSystem::getNEqs() const  {
    TRACE(14,"TaSystem::getNDofs()");
    us Neqs=0;
    for(us i=0;i<getNSegs();i++)
      Neqs+=segs.at(i)->getNEqs();
    return Neqs;
  }
  void TaSystem::connectSegs(us seg1,us seg2,SegCoupling sc){
    TRACE(14,"TaSystem::ConnectSegs()");
    // Basic check if nothing is wrong
    if(max(seg1,seg2)>=getNSegs())
      {
	WARN("Segment number is higher than available number of segments. ");
	return;
      }
    segConnections.push_back(SegConnection(seg1,seg2,sc));
    hasInit=false;
  }

  void TaSystem::show(us showvertices){
    checkInit();
    TRACE(14,"checkInit() done");
    cout << "########################## Showing TaSystem...\n";
    cout << "Showing Global configuration...\n";
    gc.show();
    cout << "Now showing segments int TaSystem...\n";
    for(us i=0;i<getNSegs();i++){
      TRACE(13,"Showing segment "<<i <<"..");
      segs[i]->show(showvertices);
    }
  }
  void TaSystem::setGc(const Globalconf& gc){
    TRACE(14,"TaSystem::setGc()");
    this->gc=gc;
    hasInit=false;
  }
  void TaSystem::jacTriplets(TripletList& trips){
    TRACE(14,"TaSystem::jacTriplets()");
    
    Jacobian jnew;
    us Nsegs=getNSegs();
    for(us j=0;j<getNSegs();j++){
      TRACE(14,"System loop, segment " << j);
      segment::SegBase& curseg=*segs[j].get();
      curseg.jac(jnew);
      TRACE(10,"Creation of Jacobian for segment "<< j << "done."<<endl);
    } // end for loop
    // TRACE(25,"Jac\n"<<jac);
    trips=jnew.getTriplets();

  }
  esdmat TaSystem::jac(d dummy){
    TRACE(14,"TaSystem::Jac()");
    checkInit();
    // Something interesting has to be done here later on to connect
    // the different segments in the sense that blocks of Jacobian
    // matrix parts have to be moved to the right place etc. To be
    // continued...
    const us& Ns=gc.Ns;
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
    us Nsegs=getNSegs();
    us segeqs;
    us starteq=0;
    TRACE(-1,"Nsegs:"<< Nsegs);
    for(us i=0;i<Nsegs;i++){
      segeqs=segs[i]->getNEqs();
      Error.subvec(starteq,starteq+segeqs-1)=segs[i]->error();
      starteq+=segeqs;
    }
    return error;
  }
  evd TaSystem::getRes(){
    TRACE(14,"TaSystem::getRes()");
    checkInit();
    us Ndofs=getNDofs();
    TRACE(14,"TaSystem::GetRes(), Ndofs:"<< Ndofs);
    const us& Ns=gc.Ns;
    evd res(Ndofs);
    
    vd Res(res.data(),Ndofs,false,false);

    us segdofs;
    us startdof=0;
    us Nsegs=getNSegs();
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
    us nsegs=getNSegs();
    assert(other.getNSegs()==nsegs);
    for(us i=0;i<nsegs;i++) {
      getSeg(i)->setRes(*other.getSeg(i));
    }
    
  }
  void TaSystem::resetHarmonics(){
    assert(!segs.empty());
    for(auto seg=segs.begin();seg!=segs.end();seg++)
      seg->get()->resetHarmonics();
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
      us Nsegs=getNSegs();
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
      Eigen::IOFormat CleanFmt(1,0," ",";\n","","","[","]");
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
  }
  TaSystem::~TaSystem() {
    TRACE(-5,"~TaSystem()");
    cleanup();
  }

} // namespace tasystem

