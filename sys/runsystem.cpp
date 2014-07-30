#include "system.h"
#include "segbase.h"

namespace tasystem{
  using segment::SegBase;
  dmat TAsystem::Jac(){
    TRACE(14,"TAsystem::Jac()");
    CheckInit();
    // Something interesting has to be done here later on to connect
    // the different segments in the sense that blocks of Jacobian
    // matrix parts have to be moved to the right place etc. To be
    // continued...
    TRACE(-1,"Ndofs:"<<Ndofs);
    const us& Ns=gc.Ns;
    us Nsegs=getNsegs();
    dmat jac(Ndofs,Ndofs);
    us cellblock=Neq*Ns;
    for(us j=0;j<getNsegs();j++){
      TRACE(14,"System loop, segment " << j);
      segment::Seg& curseg=*segs[j].get();
      dmat segjac=curseg.Jac();
      TRACE(14,"Obtaining sub-Jacobian for segment "<< j << "...");
      us thisndofs=curseg.getNdofs();
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
      	TRACE(14,"Coupling to left segment..");
	const SegBase& left=*(curseg.getLeft()[0]);
      	// Couple Jacobian terms
      	us othernr=left.getNumber();
      	us firstcolother=segfirstdof(othernr);
	// TRACE(100,"Firstcol other segment:"<< firstcolother);	
      	us otherndofs=segndofs(othernr);
      	// Find out if other segment is coupled to the left, or to the right
      	if(left.getRight()[0]->getNumber()==thisnr){
      	  TRACE(14,"tail of left segment "<< thisnr << " coupled to head of current segment");
      	  jac.submat(frow,firstcolother+otherndofs-cellblock,lrow,firstcolother+otherndofs-1)= \
      	    segjac.cols(0,cellblock-1);
	  // TRACE(100,"Coupling jacobian terms:\n"<<segjac.cols(0,cellblock-1));
      	}    
      	else{			// headhead coupling

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
      	  // tail of left segment coupled to head of current segment
      	  jac.submat(frow,firstcolother,lrow,firstcolother+cellblock-1)=	\
      	    segjac.cols(segjac.n_cols-cellblock,segjac.n_cols-1);
	  // TRACE(100,"Coupling jacobian terms:\n"<<segjac.cols(segjac.n_cols-cellblock,segjac.n_cols-1));
      	}    
      	else{			// headhead coupling

      	} // 
      }	  // curseg.Right()!=NULL
      TRACE(-1,"Creation of Jacobian for segment "<< j << "done."<<endl);
    } // end for loop
    return jac;
  }
  vd TAsystem::Error(){
    TRACE(14,"TAsystem::Error()");
    CheckInit();
    vd Error(Ndofs);
    us Nsegs=getNsegs();
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
  vd TAsystem::GetRes(){
    TRACE(14,"TAsystem::GetRes(), Ndofs:"<< Ndofs);
    CheckInit();
    const us& Ns=gc.Ns;
    vd Res(Ndofs);
    us segdofs;
    us startdof=0;
    us Nsegs=getNsegs();
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
    if(Res.size()==Ndofs){
      us segdofs;
      us startdof=0;
      us Nsegs=getNsegs();
      for(us i=0;i<Nsegs;i++){
	segdofs=segs[i]->getNdofs();
	segs[i]->SetRes(Res.subvec(startdof,startdof+segdofs-1));
	startdof=startdof+segdofs;
      }
    }
    else
      {
	WARN("Amount of DOFS in result vector does not match system size!");
      }
  } // TAsystem::SetRes()
} //namespace tasystem
