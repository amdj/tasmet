#include "system.h"


namespace tasystem{
  inline us max(us s,us t){  return s? s>=t : t;}

  TAsystem::TAsystem(const Globalconf& gc):gc(gc){
    TRACE(14,"TAsystem::TAsystem(gc)");
  }
  TAsystem::TAsystem(const TAsystem& o)
  {
    TRACE(14,"TAsystem::TAsystem(TAsystem&)");
    copyTAsystem(o);
  }
  void TAsystem::copyTAsystem(const TAsystem& o){
    TRACE(14,"TAsystem::copyTAsystem()");
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
  TAsystem& TAsystem::operator=(const TAsystem& other){
    TRACE(14,"TAsystem::operator=()");
    cleanup();
    copyTAsystem(other);
    return *this;
  }
  void TAsystem::cleanup(){
    segs.clear();
    segConnections.clear();
    segfirstdof.zeros();
    segndofs.zeros();
    hasInit=false;
  }
  void TAsystem::addSeg(const SegBase& seg){
    TRACE(14,"TAsystem::addseg()");
    hasInit=false;
    segs.emplace_back(seg.copy());
    segs[getNSegs()-1]->setNumber(getNSegs()-1);
  }
  SegBase* TAsystem::getSeg(us i) const { return (*this)[i];}
  
  SegBase* TAsystem::operator[](us i) const {
    us nSegs=getNSegs();
    if(i<nSegs)
      return segs[i].get();
    else
      return NULL;
  }

  void TAsystem::init(){
    TRACE(14,"TAsystem::init()");
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
  TAsystem::~TAsystem() {
    TRACE(-5,"~TAsystem()");
    cleanup();
  }

  
  us TAsystem::getNDofs()
  {
    TRACE(14,"TAsystem::getNDofs()");
    us Ndofs=0;
    for(us i=0;i<getNSegs();i++)
      Ndofs+=segs.at(i)->getNDofs();
    return Ndofs;
  }
  void TAsystem::connectSegs(us seg1,us seg2,SegCoupling sc){
    TRACE(14,"TAsystem::ConnectSegs()");
    // Basic check if nothing is wrong
    if(max(seg1,seg2)>=getNSegs())
      {
	WARN("Segment number is higher than available number of segments. ");
	return;
      }
    segConnections.push_back(SegConnection(seg1,seg2,sc));
    hasInit=false;
  }

  void TAsystem::show(bool showvertices){
    cout << "########################## Showing TAsystem...\n"		\
      ;
    checkInit();
    TRACE(14,"checkInit() done");
    gc.show();
    for(us i=0;i<getNSegs();i++){
      TRACE(13,"Showing segment "<<i <<"..");
      segs[i]->show(showvertices);
    }
  }
  void TAsystem::checkInit(){
    TRACE(14,"TAsystem::CheckInit()");
    if(!hasInit){
      init();
      hasInit=true;
    }
  }
  void TAsystem::setGc(const Globalconf& gc){
    TRACE(14,"TAsystem::setGc()");
    this->gc=gc;
    hasInit=false;
  }
  // void TAsystem::setnodes(us segnr,us nl,us nr){
  //   TRACE(14,"TAsystem::setnodes");
  //   assert(segnr>0 && segnr<Nsegs);
  //   segs[segnr]->setnodes(nl,nr);
  // }

} // namespace tasystem

