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

    dmat jac(Ndofs,Ndofs,fillwith::zeros);
    for(us j=0;j<Nsegs;j++){
      dmat segjac=segs[j]->Jac();
      TRACE(0,"Seg jacobian size:"<<segjac.n_cols);
      us thisndofs=gc.Ns*Neq*segs[j]->Ncells;
      us frow=j*thisndofs;
      us fcol=j*thisndofs;
      us lrow=(j+1)*thisndofs-1;
      us lcol=(j+1)*thisndofs-1;
      jac.submat(frow,fcol,lrow,lcol)=segjac;
    }
    return jac;
  }
  


  void TAsystem::addseg(Seg& s){
    TRACE(0,"TAsystem::addseg()");
    segs.push_back(&s);
    if(Nsegs==0){
      startdof.push_back(0);
      enddof.push_back(s.Ncells*Neq*gc.Ns-1);
    }
    else{
      WARN("addseg implementation incomplete!");
    }
    Nsegs++;
    Ndofs+=s.Ncells*gc.Ns*Neq;
  }
  Seg& TAsystem::operator[](us i){
    return *(segs.at(i));
  }
  TAsystem::~TAsystem(){}
} // namespace tasystem
