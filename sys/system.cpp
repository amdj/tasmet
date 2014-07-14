#include "system.h"
#include "impedancebc.h"



namespace tasystem{

  void copyallsegsbc(TAsystem& to,const TAsystem& from){
    TRACE(14,"copyallsegsbc()");
    for(us i=0;i<from.getNsegs();i++)
      if(from[i]!=NULL)
	to.addseg(*from[i]);
      else
	TRACE(14,"Warning! Segment is NULL");
    for(us i=0;i<from.getNbc();i++)
      if(from.getBc(i)!=NULL)
	to.addbc(*from.getBc(i));
      else
	TRACE(14,"Warning! Bc is NULL");
    
  }
  us computeNdofs(const TAsystem& t)
  {
    us Ndofs=0; 
    for(us i=0;i<t.getNsegs();i++)
      Ndofs+=t[i]->getNdofs();
    return Ndofs;
  }
  // void copybc()

  TAsystem::TAsystem(const Globalconf& gc):gc(gc){}
  TAsystem::TAsystem(const TAsystem& o): gc(o.gc)
  {
    TRACE(14,"TAsystem copy cc");
    copyallsegsbc(*this,o);
  }

  TAsystem& TAsystem::operator=(const TAsystem& other){
    cleanup();
    gc=other.gc;
    copyallsegsbc(*this,other);
    return *this;
  }
  void TAsystem::show(){
    cout << "Showing TAsystem...\n"		\
      ;
    gc.show();
    for(us i=0;i<Nsegs;i++)
      segs[i]->show();
  }
  void TAsystem::cleanup(){
    for(us i=0;i<Nsegs;i++)
      delete segs[i];
    for(us i=0;i<Nbc;i++)
      delete bcvertices[i];
    Nsegs=0;
    Ndofs=0;
    Nbc=0;
  }
  BcVertex* TAsystem::getBc(us i) const {
      if(i<Nbc)
	return bcvertices[i];
      else
	return NULL;
	
    }
  void TAsystem::addseg(const Seg& seg){
    TRACE(14,"TAsystem::addseg()");
    segs.push_back(copyseg(seg));
    Nsegs++;			// Update number of segments
    Ndofs=computeNdofs(*this);
  }
  void TAsystem::addbc(const BcVertex& vertex){
    TRACE(14,"TAsystem::addbc()");
    bcvertices.push_back(copybc(vertex));
    Nbc++;
  }
  void TAsystem::CheckInit(){
    if(!hasinit){
      Init();
      hasinit=true;
    }
  }
  void TAsystem::Init(){
    TRACE(14,"TAsystem::Init()");
    for(us i=0;i<Nbc;i++)
      {
	TRACE(9,"Connecting boundary condition "<<i<<"...");
	connectbc(*segs[bcvertices[i]->segNumber()],*bcvertices[i]);
      }

    for(us i=0;i<Nsegs;i++)
      {
	TRACE(9,"Initializing Segment "<<i<<"...");
	segs[i]->Init(gc);
      }
    Ndofs=computeNdofs(*this);
  }
  void TAsystem::setGc(const Globalconf& gc){
    this->gc=gc;
    for(us i=0;i<Nsegs;i++){
      segs[i]->Init(gc);
    }
  }
  void TAsystem::setnodes(us segnr,us nl,us nr){
    TRACE(14,"TAsystem::setnodes");
    assert(segnr>0 && segnr<Nsegs);
    segs[segnr]->setnodes(nl,nr);
  }
  vd TAsystem::GetRes(){
    TRACE(14,"TAsystem::GetRes(), Ndofs:"<< Ndofs);
    CheckInit();
    const us& Ns=gc.Ns;
    vd Res(Ndofs);
    us segdofs;
    us startdof=0;
    for(us i=0;i<Nsegs;i++){
      segdofs=segs[i]->getNdofs();
      Res.subvec(startdof,startdof+segdofs-1)=segs[i]->GetRes();
      startdof=startdof+segdofs;
      TRACE(4,"Seg:"<<i<<", Ndofs: "<<segdofs);

    }
    return Res;
  }
  void TAsystem::SetRes(vd Res){
    CheckInit();
    TRACE(14,"TAsystem::SetRes(vd res)");
    us segdofs;
    us startdof=0;

    for(us i=0;i<Nsegs;i++){
      segdofs=segs[i]->getNdofs();
      segs[i]->SetRes(Res.subvec(startdof,startdof+segdofs-1));
      startdof=startdof+segdofs;
    }
  }
  dmat TAsystem::Jac(){
    TRACE(14,"TAsystem::operator()() return Jacobian matrix");
    CheckInit();
    // Something interesting has to be done here later on to connect
    // the different segments in the sense that blocks of Jacobian
    // matrix parts have to be moved to the right place etc. To be
    // continued...
    TRACE(-1,"Ndofs:"<<Ndofs);
    const us& Ns=gc.Ns;
    dmat jac(Ndofs,Ndofs);
    us cellblock=Neq*Ns;
    for(us j=0;j<Nsegs;j++){
      TRACE(14,"back in system loop, j=" << j);
      segment::Seg& curseg=*segs[j];
      dmat segjac=curseg.Jac();
      TRACE(14,"back in system loop");
      us thisndofs=curseg.getNdofs();
      us frow=j*thisndofs;
      us fcol=j*thisndofs;	 // First col
      us lrow=(j+1)*thisndofs-1; // last row
      us lcol=(j+1)*thisndofs-1;
      TRACE(14,"Filling system Jacobian submat...");
      // cout << "Segjac:\n"<<segjac;
      jac.submat(frow,fcol,lrow,lcol)=			\
	segjac.cols(Neq*gc.Ns,segjac.n_cols-1-Neq*Ns);

      TRACE(14,"Jacobian submat succesfully filled.");

      // if(curseg.Left()[0]!=NULL){
      // 	TRACE(14,"Coupling to left segment..");
      // 	// Couple Jacobian terms
      // 	us othernr=curseg.Left()[0]->getNumber();
      // 	us firstcol=segfirstcol(othernr);
      // 	us otherndofs=segndofs(othernr);
      // 	// Find out if other segment is coupled to the left, or to the right
      // 	if(*(curseg.Left()[0]->Right()[0])==curseg){
      // 	  // tail of left segment coupled to head of current segment
      // 	  jac.submat(frow,firstcol+otherndofs-cellblock,lrow,firstcol+otherndofs-1)= \
      // 	    segjac.cols(0,cellblock-1);
      // 	}    
      // 	else{			// headhead coupling

      // 	} // 
      // }	  // curseg.Left()!=NULL

      
      // if(curseg.Right()[0]!=NULL){
      // 	// Couple Jacobian terms
      // 	TRACE(14,"Coupling to right segment..");
      // 	us othernr=curseg.Right()[0]->getNumber();
      // 	us firstcol=segfirstcol(othernr);
      // 	us otherndofs=segndofs(othernr);
      // 	// Find out if other segment is coupled to the left, or to the right
      // 	if(*(curseg.Right()[0]->Left()[0])==curseg){
      // 	  // tail of left segment coupled to head of current segment
      // 	  jac.submat(frow,firstcol,lrow,firstcol+cellblock-1)=	\
      // 	    segjac.cols(segjac.n_cols-cellblock,segjac.n_cols-1);
      // 	}    
      // 	else{			// headhead coupling

      // 	} // 
      // }	  // curseg.Right()!=NULL
      TRACE(-1,"Creation of Jacobian for segment "<< j << "done."<<endl);
    } // end for loop
    return jac;
  }
  vd TAsystem::Error(){
    TRACE(14,"TAsystem::Error()");
    CheckInit();
    vd Error(Ndofs);
    us segdofs;
    us startdof=0;

    TRACE(-1,"Nsegs:"<< Nsegs);
    for(us i=0;i<Nsegs;i++){
      segdofs=segs[i]->getNdofs();
      Error.subvec(startdof,startdof+segdofs-1)=segs[i]->Error();
      startdof=startdof+segdofs;
    }
    return Error;
  }

  Seg* TAsystem::operator[](us i) const {
    if(i<Nsegs)
      return (segs[i]);
    else
      return NULL;
  }
  
  TAsystem::~TAsystem() {
    TRACE(-5,"~TAsystem()");
    cleanup();
  }
} // namespace tasystem

