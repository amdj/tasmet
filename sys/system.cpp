#include "system.h"

namespace tasystem{

  TAsystem::TAsystem(Globalconf& g):gc(g),Ns(g.Ns){
    TRACE(1,"TAsystem::TAsystem(SegArray)");
    Nsegs=0;
    Ndofs=0;
    
  }
  vd TAsystem::Error(){
    TRACE(1,"TAsystem::Error()");
    
    vd Error(Ndofs);
    us segdofs;
    us startdof=0;

    for(us i=0;i<Nsegs;i++){
      segdofs=segndofs(i);
      startdof=segfirstcol(i);
      Error.subvec(startdof,startdof+segdofs-1)=segs[i]->Error();
      startdof=startdof+segdofs;
    }
    return Error;
  }
  void TAsystem::setnodes(us segnr,us nl,us nr){
    TRACE(1,"TAsystem::setnodes");
    assert(segnr>0 && segnr<Nsegs);
    segs[segnr]->setnodes(nl,nr);
  }
  vd TAsystem::GetRes(){
    TRACE(1,"TAsystem::GetRes()");
    vd Res(Ndofs);
    us segdofs;
    us startdof=0;
    for(us i=0;i<Nsegs;i++){
      segdofs=segs[i]->getNcells()*Ns*Neq;
      Res.subvec(startdof,startdof+segdofs-1)=segs[i]->GetRes();
      startdof=startdof+segdofs;
    }
    return Res;
  }
  void TAsystem::SetRes(vd Res){
    TRACE(1,"TAsystem::SetRes(vd res)");
    us segdofs;
    us startdof=0;

    for(us i=0;i<Nsegs;i++){
      segdofs=segs[i]->getNcells()*Ns*Neq;
      segs[i]->SetRes(Res.subvec(startdof,startdof+segdofs-1));
    }
  }
  dmat TAsystem::Jac(){
    TRACE(1,"TAsystem::operator()() return Jacobian matrix");
    // Something interesting has to be done here later on to connect
    // the different segments in the sense that blocks of Jacobian
    // matrix parts have to be moved to the right place etc. To be
    // continued...
    TRACE(-1,"Ndofs:"<<Ndofs);
    const us& Ns=gc.Ns;
    dmat jac(Ndofs,Ndofs);
    us cellblock=Neq*Ns;
    for(us j=0;j<Nsegs;j++){
      segment::Seg& curseg=*segs[j];
      dmat segjac=curseg.Jac();
      us thisndofs=gc.Ns*Neq*curseg.getNcells();
      us frow=j*thisndofs;
      us fcol=j*thisndofs;	 // First col
      us lrow=(j+1)*thisndofs-1; // last row
      us lcol=(j+1)*thisndofs-1;
      jac.submat(frow,fcol,lrow,lcol)=		\
	segjac.cols(Neq*gc.Ns,segjac.n_cols-Neq*Ns);

      if(curseg.Left()[0]!=NULL){
	TRACE(10,"Coupling to left segment..");
	// Couple Jacobian terms
	us othernr=curseg.Left()[0]->getNumber();
	us firstcol=segfirstcol(othernr);
	us otherndofs=segndofs(othernr);
	// Find out if other segment is coupled to the left, or to the right
	if(*(curseg.Left()[0]->Right()[0])==curseg){
	  // tail of left segment coupled to head of current segment
	  jac.submat(frow,firstcol+otherndofs-cellblock,lrow,firstcol+otherndofs-1)= \
	    segjac.cols(0,cellblock-1);
	}    
	else{			// headhead coupling

	} // 
      }	  // curseg.Left()!=NULL

      
      if(curseg.Right()[0]!=NULL){
	// Couple Jacobian terms
	TRACE(10,"Coupling to right segment..");
	us othernr=curseg.Right()[0]->getNumber();
	us firstcol=segfirstcol(othernr);
	us otherndofs=segndofs(othernr);
	// Find out if other segment is coupled to the left, or to the right
	if(*(curseg.Right()[0]->Left()[0])==curseg){
	  // tail of left segment coupled to head of current segment
	  jac.submat(frow,firstcol,lrow,firstcol+cellblock-1)=	\
	    segjac.cols(segjac.n_cols-cellblock,segjac.n_cols-1);
	}    
	else{			// headhead coupling

	} // 
      }	  // curseg.Right()!=NULL
      
    } // end for loop
    return jac;
  }

  void TAsystem::addseg(Seg& s){
    TRACE(10,"TAsystem::addseg()");
    if(Nsegs==MAXSEGS){
      TRACE(0,"Warning: maximum segments reached");
      return;
    }
    // Put the segment in the array
    segs.push_back(&s);
    s.Init();
    segfirstcol(Nsegs)=Ndofs;
    Ndofs+=s.getNcells()*gc.Ns*Neq;	// number of cells times number of equations times number of time samples
    segndofs(Nsegs)=(s.getNcells()*gc.Ns*Neq);
    Nsegs++;			// Update number of segments    
  }
  Seg& TAsystem::operator[](us i){
    return *(segs.at(i));
  }

} // namespace tasystem

