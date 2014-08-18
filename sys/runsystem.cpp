#include "system.h"


namespace tasystem{
  using segment::SegBase;

  dmat taSystem::jac(){
    TRACE(14,"taSystem::Jac()");
    checkInit();
    // Something interesting has to be done here later on to connect
    // the different segments in the sense that blocks of Jacobian
    // matrix parts have to be moved to the right place etc. To be
    // continued...
    const us& Ns=gc.Ns;
    us Ndofs=getNDofs();
    us Nsegs=getNSegs();
    TRACE(-1,"Ndofs:"<<Ndofs);
    dmat jac(Ndofs,Ndofs,fillwith::zeros);
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
    // TRACE(15,"Jac\n"<<jac);
    return jac;
  }
  vd taSystem::error(){
    TRACE(14,"taSystem::Error()");
    checkInit();
    us Ndofs=getNDofs();
    vd Error(getNDofs());
    us Nsegs=getNSegs();
    us segdofs;
    us startdof=0;
    TRACE(-1,"Nsegs:"<< Nsegs);
    for(us i=0;i<Nsegs;i++){
      segdofs=segs[i]->getNDofs();
      Error.subvec(startdof,startdof+segdofs-1)=segs[i]->error();
      startdof=startdof+segdofs;
    }
    return Error;
  }
  vd taSystem::getRes(){
    TRACE(14,"taSystem::getRes()");
    checkInit();
    us Ndofs=getNDofs();
    TRACE(14,"taSystem::GetRes(), Ndofs:"<< Ndofs);
    const us& Ns=gc.Ns;
    vd Res(Ndofs);
    us segdofs;
    us startdof=0;
    us Nsegs=getNSegs();
    for(us i=0;i<Nsegs;i++){
      segdofs=segs[i]->getNDofs();
      Res.subvec(startdof,startdof+segdofs-1)=segs[i]->getRes();
      startdof=startdof+segdofs;
      TRACE(4,"Seg:"<<i<<", Ndofs: "<<segdofs);

    }
    return Res;
  }
  void taSystem::setRes(vd Res){
    checkInit();
    TRACE(14,"taSystem::SetRes(vd res)");
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
  } // taSystem::SetRes()
} //namespace tasystem
