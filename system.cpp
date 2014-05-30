#include "system.h"

namespace tasystem{

  TAsystem::TAsystem(Globalconf& g):gc(g),Ns(g.Ns){
    TRACE(0,"TAsystem::TAsystem(SegArray)");
    Nsegs=0;
    Ndofs=0;
  }
  vd TAsystem::Error(){
    TRACE(0,"TAsystem::Error()");
    
    vd Error(Ndofs);
    for(us i=0;i<Nsegs;i++){
      Error.subvec(startdof.at(i),enddof.at(i))=segs[i]->Error();
    }
    return Error;
  }
  void TAsystem::setnodes(us segnr,us nl,us nr){
    TRACE(0,"TAsystem::setnodes");
    assert(segnr>0 && segnr<Nsegs);
    segs[segnr]->setnodes(nl,nr);
  }
  vd TAsystem::GetRes(){
    TRACE(0,"TAsystem::GetRes()");
    vd Res(Ndofs);
    for(us i=0;i<Nsegs;i++){
      Res.subvec(startdof.at(i),enddof.at(i))=segs[i]->GetRes();
    }
    return Res;
  }
  void TAsystem::SetRes(vd Res){
    TRACE(0,"TAsystem::SetRes(vd res)");
    for(us i=0;i<Nsegs;i++){
      segs[i]->SetRes(Res.subvec(startdof.at(i),enddof.at(i)));
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
    TRACE(0,"TAsystem::addseg()");
    // Put the segment in the array
    segs.push_back(&s);
    Ndofs+=s.Ncells*gc.Ns*Neq;	// number of cells times number of equations times number of time samples
    if(Nsegs==0){
      startdof.push_back(0);
    }
    else{
      startdof.push_back(startdof[Nsegs]);
    }
    enddof.push_back(s.Ncells*Neq*gc.Ns-1);
    Nsegs++;			// Update number of segments    
  }
  Seg& TAsystem::operator[](us i){
    return *(segs.at(i));
  }

} // namespace tasystem

