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
      segdofs=segs[i]->Ncells*Ns*Neq;
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
      segdofs=segs[i]->Ncells*Ns*Neq;
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
      segdofs=segs[i]->Ncells*Ns*Neq;
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
    for(us j=0;j<Nsegs;j++){
      segment::Seg& curseg=*segs[j];
      dmat segjac=curseg.Jac();
      us thisndofs=gc.Ns*Neq*curseg.Ncells;
      us frow=j*thisndofs;
      us fcol=j*thisndofs;	 // First col
      us lrow=(j+1)*thisndofs-1; // last row
      us lcol=(j+1)*thisndofs-1;
      jac.submat(frow,fcol,lrow,lcol)=segjac.cols(Neq*gc.Ns,segjac.n_cols-Neq*Ns);
    }
    return jac;
  }

  void TAsystem::addseg(Seg& s){
    TRACE(1,"TAsystem::addseg()");
    // Put the segment in the array
    segs.push_back(&s);
    Ndofs+=s.Ncells*gc.Ns*Neq;	// number of cells times number of equations times number of time samples
    // if(Nsegs==0){
      // startdof.push_back(0);
    // }
    // else{
      // startdof.push_back(enddof[Nsegs]+1);
    // }
    // enddof.push_back(s.Ncells*Neq*gc.Ns-1);
    Nsegs++;			// Update number of segments    
  }
  Seg& TAsystem::operator[](us i){
    return *(segs.at(i));
  }

} // namespace tasystem

