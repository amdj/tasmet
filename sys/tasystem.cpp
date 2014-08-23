#include "tasystem.h"
#include "arma_eigen.h"


namespace tasystem{
  using segment::SegBase;
  using math_common::armaView;
  inline us max(us s,us t){  return s? s>=t : t;}

  TaSystem::TaSystem(const Globalconf& gc):gc(gc){
    TRACE(14,"TaSystem::TaSystem(gc)");
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
    segfirstdof.zeros();
    segndofs.zeros();
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
    // Last, but not leas: initialize a pointer to this tasystem in
    // globalconf
    gc.setSys(this);
    
    hasInit=true;
  }
  

  
  us TaSystem::getNDofs()
  {
    TRACE(14,"TaSystem::getNDofs()");
    us Ndofs=0;
    for(us i=0;i<getNSegs();i++)
      Ndofs+=segs.at(i)->getNDofs();
    return Ndofs;
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

  void TaSystem::show(bool showvertices){
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
  void TaSystem::checkInit(){
    TRACE(14,"TaSystem::CheckInit()");
    if(!hasInit){
      init();
      hasInit=true;
    }
  }
  void TaSystem::setGc(const Globalconf& gc){
    TRACE(14,"TaSystem::setGc()");
    this->gc=gc;
    hasInit=false;
  }

  esdmat TaSystem::jac(){
    TRACE(14,"TaSystem::Jac()");
    checkInit();
    // Something interesting has to be done here later on to connect
    // the different segments in the sense that blocks of Jacobian
    // matrix parts have to be moved to the right place etc. To be
    // continued...
    const us& Ns=gc.Ns;
    us Ndofs=getNDofs();
    us Nsegs=getNSegs();
    Eigen::MatrixXd eigjac(Ndofs,Ndofs);
    
    TRACE(-1,"Ndofs:"<<Ndofs);
    dmat jac(eigjac.data(),Ndofs,Ndofs,false);
    us cellblock=Neq*Ns;
    for(us j=0;j<getNSegs();j++){
      TRACE(14,"System loop, segment " << j);
      segment::SegBase& curseg=*segs[j].get();
      dmat&& segjac=curseg.jac();
      TRACE(14,"Obtaining sub-Jacobian for segment "<< j << "...");
      us thisndofs=curseg.getNDofs();
      us frow=segfirstdof(j);
      us fcol=segfirstdof(j);	 // First col
      TRACE(14,"First dof:"<< frow);
      us lrow=frow+thisndofs-1; // last row
      us lcol=fcol+thisndofs-1;

      // cout << "Segjac:\n"<<segjac;
      jac.submat(frow,fcol,lrow,lcol)=			\
	segjac.cols(Neq*gc.Ns,segjac.n_cols-1-Neq*Ns);

      TRACE(14,"Jacobian submat succesfully filled.");
      us thisnr=curseg.getNumber();
      if(curseg.getLeft().size()>0){
	const SegBase& left=*(curseg.getLeft()[0]);
      	// Couple Jacobian terms
      	us othernr=left.getNumber();
      	us firstcolother=segfirstdof(othernr);
	// TRACE(100,"Firstcol other segment:"<< firstcolother);	
      	us otherndofs=segndofs(othernr);
      	// Find out if other segment is coupled to the left, or to the right
      	if(left.getRight()[0]->getNumber()==thisnr){
      	  TRACE(14,"Tail of left segment "<< othernr << " coupled to head of segment " << thisnr<<".");
      	  jac.submat(frow,firstcolother+otherndofs-cellblock,lrow,firstcolother+otherndofs-1)= \
      	    segjac.cols(0,cellblock-1);
      	}    
      	else{			// headhead coupling
	  WARN("Head-head coupling not yet implemented! Exiting.");
	  exit(1);
      	} // 
      }	  // curseg.Left()!=NULL

      
      if(curseg.getRight().size()>0){
      	// Couple Jacobian terms
	const SegBase& right=*curseg.getRight()[0];
      	TRACE(14,"Coupling to right segment..");
      	us othernr=right.getNumber();
      	us firstcolother=segfirstdof(othernr);
	// TRACE(100,"Firstcol other segment:"<< firstcolother);
      	us otherndofs=segndofs(othernr);
      	// Find out if other segment is coupled to the left, or to the right
      	if(curseg.getRight()[0]->getLeft()[0]->getNumber()==thisnr){
	  TRACE(14,"Tail of segment "<< thisnr << " connected to head of seg " << right.getNumber() << ".");
      	  jac.submat(frow,firstcolother,lrow,firstcolother+cellblock-1)=	\
      	    segjac.cols(segjac.n_cols-cellblock,segjac.n_cols-1);
      	}    
      	else{			// headhead coupling
	  WARN("Head-head coupling not yet implemented! Exiting.");
	  exit(1);
      	} // 
      }	  // curseg.Right()!=NULL
      TRACE(-1,"Creation of Jacobian for segment "<< j << "done."<<endl);
    } // end for loop
    // TRACE(25,"Jac\n"<<jac);
    esdmat eigsjac=eigjac.sparseView();
    eigsjac.makeCompressed();
    return eigsjac;
  }
  evd TaSystem::error(){
    TRACE(14,"TaSystem::Error()");
    checkInit();
    us Ndofs=getNDofs();
    evd error(getNDofs());
    vd Error(error.data(),getNDofs(),false,false);
    Error.zeros();
    us Nsegs=getNSegs();
    us segdofs;
    us startdof=0;
    TRACE(-1,"Nsegs:"<< Nsegs);
    for(us i=0;i<Nsegs;i++){
      segdofs=segs[i]->getNDofs();
      Error.subvec(startdof,startdof+segdofs-1)=segs[i]->error();
      startdof=startdof+segdofs;
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
      startdof=startdof+segdofs;
      TRACE(4,"Seg:"<<i<<", Ndofs: "<<segdofs);

    }
    return res;
  }
  void TaSystem::setRes(vd Res){
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
	WARN("Amount of DOFS in result vector does not match system size!");
      }
  } // TaSystem::SetRes()
  d TaSystem::getCurrentMass() const{
    TRACE(10,"TaSystem::getCurrentMass()");
    d mass=0;
    us nsegs=segs.size();
    for(us i=0;i<nsegs;i++){
      mass+=segs[i]->getCurrentMass();
    } // for loop
    return mass;
  }
  
  TaSystem::~TaSystem() {
    TRACE(-5,"~TaSystem()");
    cleanup();
  }

} // namespace tasystem

