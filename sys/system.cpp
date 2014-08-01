#include "system.h"
#include "impedancebc.h"



namespace tasystem{
  inline us max(us s,us t){  return s? s>=t : t;}

  TAsystem::TAsystem(const Globalconf& gc):gc(gc){
    TRACE(14,"TAsystem::TAsystem(gc)");
    segfirstdof.zeros();
    segndofs.zeros();
  }
  TAsystem::TAsystem(const TAsystem& o) :TAsystem(o.gc)
  {
    TRACE(14,"TAsystem::TAsystem(TAsystem&)");
    cleanup();
    TRACE(14,"TAsystem copy cc");
    copyallsegsbc(*this,o);
    segConnections=o.segConnections;
    hasInit=false;
  }
  void TAsystem::init(){
    TRACE(14,"TAsystem::init()");
    us Nsegs=getNSegs();
    us Nbc=getNBc();
    for(us i=0;i<Nsegs;i++)
      {
	TRACE(9,"Initializing Segment "<<i<<"...");
	segs[i]->init(gc);
      }
    for(us i=0;i<Nbc;i++)
      {
    	TRACE(9,"Connecting boundary condition "<<i<<"...");
    	connectbc(*segs[bcvertices[i]->segNumber()],*bcvertices[i]);
      }
    us i=0;
    for(auto v=segConnections.begin();v!=segConnections.end();++v){
      TRACE(90,"Connecting segment connection " << i << "...");
      coupleSegs(*v,*this);
      i++;
    }

    // Now the weight functions are not yet updated for the vertices
    // that have to be connected. For now we see no other option than
    // running seg->init again
    for(us i=0;i<Nsegs;i++)
      {
	TRACE(9,"Initializing Segment "<<i<<"...");
	segs[i]->init(gc);
      }
    
    
    computeNdofs();
    if(Ndofs>MAXNDOFS)
      {
	WARN("Way too many DOFS required: Ndofs=" <<Ndofs << ". Exiting...\n");
	exit(1);
      }
    hasInit=true;
  }

  TAsystem& TAsystem::operator=(const TAsystem& other){
    cleanup();
    gc=other.gc;
    copyallsegsbc(*this,other);
    segConnections=other.segConnections;
    hasInit=false;
    return *this;
  }
  TAsystem::~TAsystem() {
    TRACE(-5,"~TAsystem()");
    cleanup();
  }

  
  void TAsystem::computeNdofs()
  {
    TRACE(14,"TAsystem::computeNdofs()");
    Ndofs=0;
    segfirstdof(0)=0;

    for(us i=0;i<getNSegs();i++)
      {
	assert(segs[i].get()!=NULL);
	us thisndofs=segs[i]->getNDofs();
	TRACE(12,"This segment ndofs:"<< thisndofs);
	Ndofs+=thisndofs;
	segndofs(i)=thisndofs;
	if(i>0) segfirstdof(i)=segfirstdof(i-1)+segndofs(i-1);
      }
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
  void TAsystem::cleanup(){
    Ndofs=0;
    segs.clear();
    bcvertices.clear();
    segConnections.clear();
    hasInit=false;
  }
  BcVertex* TAsystem::getBc(us i) const {
    us Nbc=getNBc();
    if(i<Nbc)
      return bcvertices[i].get();
    else
      return NULL;
	
  }
  void TAsystem::addSeg(const Seg& seg){
    TRACE(14,"TAsystem::addseg()");
    hasInit=false;
    segs.emplace_back((copyseg(seg)));
    segs[getNSegs()-1]->setNumber(getNSegs()-1);
  }
  void TAsystem::addBc(const BcVertex& vertex){
    TRACE(14,"TAsystem::addbc()");
    bcvertices.emplace_back(copybc(vertex));
    hasInit=false;
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

  Seg* TAsystem::operator[](us i) const {
    us nSegs=getNSegs();
    
    if(i<nSegs)
      return (segs[i].get());
    else
      return NULL;
  }
  
} // namespace tasystem

