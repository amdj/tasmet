#include "seg.h"

namespace segment{
   

    
  // Vertex* vertexfrombc(BcVertex* orig){
  //   if(orig->gettype().compare("RightImpedance")==0)
  //     return static_cast<Vertex*>(static_cast<RightImpedance*>(orig));
  //   else if(orig->gettype().compare("TwImpedance")==0)
  //     return static_cast<Vertex*>(static_cast<TwImpedance*>(orig));
  //   else if(orig->gettype().compare("RightIsoTWall")==0)
  //     return static_cast<Vertex*>(static_cast<RightIsoTWall*>(orig));
  //   else if(orig->gettype().compare("LeftPressure")==0)
  //     return static_cast<Vertex*>(static_cast<LeftPressure*>(orig));
  //   else{
  //     WARN("Boundary condition type not understood");
  //     abort();
  //   }
  // }
  
 
  Seg::Seg(const Geom& geom):SegBase(geom){
    TRACE(13,"Seg::Seg(Geom)");
    // The Jacobian matrix is larger than the number of dofs for the
    // connection terms other segments
    // us& Ns=gc.Ns;
  }
  Seg::Seg(const Seg& other): SegBase(other){}
  Seg& Seg::operator=(const Seg& o){
    SegBase::operator=(o);
    return *this;
  }
  void Seg::init(const tasystem::Globalconf& gc){
    TRACE(13,"Seg::init()");
    SegBase::init(gc);
  } // Seg::Init
  void Seg::show(us showVertices) const {
    TRACE(18,"Seg::show()");
    cout << "Showing segment of type " << getType() <<" with number "<<getNumber()<< ".\n";
    cout << "Number: "<< getNumber() << ".\n";
    cout << "Geometry: \n";
    geom.show();
    assert(vvertex.size()!=0);
    for(auto s=getLeft().begin();s!=getLeft().end();s++)
	cout << "Left segment:" << *s << "\n";
    for(auto s=getRight().begin();s!=getRight().end();s++)
	cout << "Right segment:" << *s << "\n";
    if(showVertices==1)
      this->showVertices();
  }
  void Seg::showVertices() const {
    for(us i=0;i<vvertex.size();i++)
      vvertex[i]->show();
  }
  void Seg::jac(Jacobian& tofill) const{			// Return Jacobian matrix of error operator
    // sdmat Seg::Jac(){			// Return Jacobian matrix of error operator    
    TRACE(8," Seg::Jac() for Segment "<< getNumber() << ".");
    const us& Ns=gc->Ns;
    us nVertex=vvertex.size();

    for(us j=0;j<nVertex;j++){			   // Fill the Jacobian
      TRACE(3,"Obtaining vertex Jacobian...");
      vvertex.at(j)->jac(tofill);
    }	// end for
    // cout <<"Segment" << getNumber() <<" Jacobian done. Jac is:\n"<< Jacobian;
    // cout << "Number of colums in this jacobian" << Jacobian.n_cols<<"\n";
    TRACE(8,"Segment Jacobian done.");
  }
  vd Seg::getRes() const {
    TRACE(8,"Seg::GetRes()");
    assert(vvertex.size()!=0);
    assert(gc!=NULL);
    vd Result(getNDofs(),fillwith::zeros);
    us nVertex=vvertex.size();    
    us Ns=gc->Ns;
    us vndofs,curpos=0;
    for(us k=0; k<nVertex;k++) {
      vndofs=vvertex.at(k)->getNDofs();
      Result.subvec(curpos,curpos+vndofs-1)=vvertex[k]->getRes();
      curpos+=vndofs;
    }
    return Result;
  }

  vd Seg::error() const{
    TRACE(8,"Seg::Error()");
    assert(vvertex.size()!=0);
    assert(gc!=NULL);
    vd error(getNEqs(),fillwith::zeros);
    us nVertex=vvertex.size();    
    us Ns=gc->Ns;
    us vneqs,curpos=0;

    for(us k=0; k<nVertex;k++) {
      vneqs=vvertex.at(k)->getNEqs();
      error.subvec(curpos,curpos+vneqs-1)=vvertex[k]->error();
      curpos+=vneqs;
    }
    return error;
  }
  vd Seg::domg() const{
    TRACE(8,"Seg::Error()");
    const us& Ns=gc->Ns;
    us nVertex=vvertex.size();    
    vd domg(getNDofs(),fillwith::zeros);
    us vndofs,curpos=0;
    for(us k=0; k<nVertex;k++) {
      vndofs=vvertex.at(k)->getNDofs();
      domg.subvec(curpos,curpos+vndofs-1)=vvertex[k]->domg();
      curpos+=vndofs;
    }
    return domg;
  }

  void Seg::setRes(vd res){
    TRACE(8,"Seg::SetRes()");
    assert(res.size()==getNDofs());
    // const us& Neq=(vvertex[0]).Neq;
    const us& Ns=gc->Ns;
    us vertexdofs;
    us firstdof=0;
    for(us k=0; k<vvertex.size();k++) {
      vertexdofs=vvertex.at(k)->getNDofs();
      vvertex.at(k)->setRes(res.subvec(firstdof,firstdof+vertexdofs-1));
      firstdof+=vertexdofs;
    }
  }
  
}		 // Namespace segment
