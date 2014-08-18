#include "system.h"


namespace tasystem{
  inline us max(us s,us t){  return s? s>=t : t;}

  taSystem::taSystem(const Globalconf& gc):gc(gc){
    TRACE(14,"taSystem::taSystem(gc)");
  }
  taSystem::taSystem(const taSystem& o)
  {
    TRACE(14,"taSystem::taSystem(taSystem&)");
    copytaSystem(o);
  }
  void taSystem::copytaSystem(const taSystem& o){
    TRACE(14,"taSystem::copytaSystem()");
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
  taSystem& taSystem::operator=(const taSystem& other){
    TRACE(14,"taSystem::operator=()");
    cleanup();
    copytaSystem(other);
    return *this;
  }
  void taSystem::cleanup(){
    segs.clear();
    segConnections.clear();
    segfirstdof.zeros();
    segndofs.zeros();
    hasInit=false;
  }
  void taSystem::addSeg(const SegBase& seg){
    TRACE(14,"taSystem::addseg()");
    hasInit=false;
    segs.emplace_back(seg.copy());
    segs[getNSegs()-1]->setNumber(getNSegs()-1);
  }
  SegBase* taSystem::getSeg(us i) const { return (*this)[i];}
  
  SegBase* taSystem::operator[](us i) const {
    us nSegs=getNSegs();
    if(i<nSegs)
      return segs[i].get();
    else
      return NULL;
  }

  void taSystem::init(){
    TRACE(14,"taSystem::init()");
    us Nsegs=getNSegs();
    us Ndofs=0;

    us i=0;
    for(auto v=segConnections.begin();v!=segConnections.end();++v){
      TRACE(90,"Connecting segment connection " << i << "...");
      coupleSegs(*v,*this);
      i++;
    }
    segfirstdof(0)=0;
    for(us i=0;i<Nsegs;i++)
      {
	TRACE(9,"Initializing Segment "<<i<<"...");
	assert(segs.at(i));
	segs.at(i)->init(gc);
	us thisndofs=segs.at(i)->getNDofs();
	TRACE(12,"Ndofs for segment "<< i << ": "<<thisndofs);

	segndofs(i)=thisndofs;
	Ndofs+=thisndofs;
	if(i>0)
	  segfirstdof(i)=segfirstdof(i-1)+segndofs(i-1);
	TRACE(12,"This segment ndofs:"<< thisndofs);
      }
    TRACE(10,"Segment initialization done. Total NDofs:"<< Ndofs);
    if(Ndofs>MAXNDOFS)
      {
    	WARN("Way too many DOFS required: Ndofs=" <<Ndofs << ". Exiting...\n");
    	exit(1);
      }
    hasInit=true;
  }
  taSystem::~taSystem() {
    TRACE(-5,"~taSystem()");
    cleanup();
  }

  
  us taSystem::getNDofs()
  {
    TRACE(14,"taSystem::getNDofs()");
    us Ndofs=0;
    for(us i=0;i<getNSegs();i++)
      Ndofs+=segs.at(i)->getNDofs();
    return Ndofs;
  }
  void taSystem::connectSegs(us seg1,us seg2,SegCoupling sc){
    TRACE(14,"taSystem::ConnectSegs()");
    // Basic check if nothing is wrong
    if(max(seg1,seg2)>=getNSegs())
      {
	WARN("Segment number is higher than available number of segments. ");
	return;
      }
    segConnections.push_back(SegConnection(seg1,seg2,sc));
    hasInit=false;
  }

  void taSystem::show(bool showvertices){
    cout << "########################## Showing taSystem...\n"		\
      ;
    checkInit();
    TRACE(14,"checkInit() done");
    gc.show();
    for(us i=0;i<getNSegs();i++){
      TRACE(13,"Showing segment "<<i <<"..");
      segs[i]->show(showvertices);
    }
  }
  void taSystem::checkInit(){
    TRACE(14,"taSystem::CheckInit()");
    if(!hasInit){
      init();
      hasInit=true;
    }
  }
  void taSystem::setGc(const Globalconf& gc){
    TRACE(14,"taSystem::setGc()");
    this->gc=gc;
    hasInit=false;
  }
  // void taSystem::setnodes(us segnr,us nl,us nr){
  //   TRACE(14,"taSystem::setnodes");
  //   assert(segnr>0 && segnr<Nsegs);
  //   segs[segnr]->setnodes(nl,nr);
  // }

} // namespace tasystem

