#include "seg.h"
#include "bcvertex.h"
#define Neq (5)

namespace segment{
  
  
  Seg::Seg(Geom geom):SegBase(geom){
    TRACE(13,"Seg::Seg(Geom)");
    nDofs=0;
    type="Seg";
    // The Jacobian matrix is larger than the number of dofs for the
    // connection terms other segments
    // us& Ns=gc.Ns;
  }
  Seg::Seg(const Seg& other): SegBase(other){
    TRACE(13,"Seg copy constructor");
    this->gc=other.gc;
  }
  void Seg::init(const tasystem::Globalconf& gc){
    TRACE(13,"Seg::init()");
    // NO Do not clear segments! Boundary conditions have been added
    // vvertex.clear();
    SegBase::init(gc);
    
    us nVertex=geom.nCells;
    nDofs=geom.nCells*gc.Ns*Neq;
    this->gc=&gc;

  } // Seg::Init
  void Seg::cleanup(){
    nDofs=0;
  }
  void Seg::show(bool showVertices){
    TRACE(13,"Seg::show()");
    geom.show();
    for(auto s=getLeft().begin();s!=getLeft().end();s++)
	cout << "Left segment:" << *s << "\n";
    for(auto s=getRight().begin();s!=getRight().end();s++)
	cout << "Right segment:" << *s << "\n";

    if(showVertices==true)
      this->showVertices();
  }
  void Seg::showVertices(){
    for(us i=0;i<vvertex.size();i++)
      vvertex[i]->show();
  }

  void Seg::setLeftBc(Vertex* v){ // The segment owns the bc from then on!
    TRACE(13,"Seg::setLeftbc()");
    us nVertex=vvertex.size();
    
    assert(nVertex>0);
    vvertex[0].reset(v);
  }

  void Seg::setRightBc(Vertex* v){
    TRACE(13,"Seg::setRightbc()");
    us nVertex=vvertex.size();
    assert(nVertex>0);
    vvertex[nVertex-1].reset(v);
  }    
  
  dmat Seg::jac() const{			// Return Jacobian matrix of error operator
    // sdmat Seg::Jac(){			// Return Jacobian matrix of error operator    
    TRACE(8," Seg::Jac() for Segment "<< getNumber() << ".");
    const us& Ns=gc->Ns;
    us nVertex=vvertex.size();
    dmat Jacobian (nVertex*Neq*Ns,(nVertex+2)*Neq*Ns,fillwith::zeros);    
    // if(Jacobian.size()==0)
      // Jacobian=dmat(nVertex*Neq*Ns,(nVertex+2)*Neq*Ns);    
    dmat vJac(Neq*Ns,3*Neq*Ns,fillwith::zeros);
    us firstrow,firstcol,lastrow,lastcol;
    TRACE(8,"Filling Segment Jacobian matrix for segment "<< getNumber() <<"...");
    // #pragma omp parallel for
    for(us j=0;j<nVertex;j++){			   // Fill the Jacobian
      TRACE(3,"Obtaining vertex Jacobian...");
      vJac=vvertex[j]->jac();
      // The row height of a vertex jacobian matrix is Neq*Ns, The
      // column with is 3*Neq*Ns, since the neigbouring vertex has to
      // be found
      firstrow=j*Neq*Ns;
      lastrow=(j+1)*Neq*Ns-1;
      TRACE(3,"j:"<<j);
      if(j>0 && j<nVertex-1){
	TRACE(0,"interior vertex jacobian, j="<<j);
	firstcol=(j)*Neq*Ns;
	lastcol=(j+3)*Neq*Ns-1;
      }
      else if(j==0){
	if(getLeft().size()==0){	// First node
	  TRACE(10,"First node not connected to other segments")
	  firstcol=Neq*Ns;
	  lastcol=4*Neq*Ns-1;
	} else{
	  TRACE(10,"First node IS connected to other segment")
	  firstcol=0;	  // Fill in
	  lastcol=3*Neq*Ns-1;	  
	}
      }	// j==0
      else {			// Last vertex
	if(getRight().size()==0){
	  TRACE(10,"Last node not connected to other segments")
	  firstcol=(nVertex-2)*Neq*Ns;
	  lastcol=(nVertex+1)*Neq*Ns-1;	  
	}
	else{
	  TRACE(10,"Last node IS connected to other segment")
	  firstcol=(nVertex-1)*Neq*Ns;
	  lastcol=Jacobian.n_cols-1;
	}
      }	// j==nVertex-1
      TRACE(3,"Filling segment Jacobian with vertex subpart...");
      Jacobian.submat(firstrow,firstcol,lastrow,lastcol)=vJac;      
    }	// end for
    TRACE(8,"Segment Jacobian done.");
    return Jacobian;
  }
  vd Seg::getRes() const {
    TRACE(8,"Seg::GetRes()");
    vd Result(nDofs,fillwith::zeros);
    us Ns=gc->Ns;
    for(us k=0; k<vvertex.size();k++)
      {
	Result.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1)=vvertex[k]->getRes();
      }
    return Result;
  }

  vd Seg::error() const{
    TRACE(8,"Seg::Error()");
    const us& Ns=gc->Ns;
    us nVertex=vvertex.size();    
    vd error(nDofs,fillwith::zeros);
    for(us k=0; k<nVertex;k++)
      {
	error.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1)=vvertex[k]->error();
      }
    return error;
  }
  void Seg::setRes(vd res){
    TRACE(8,"Seg::SetRes()");
    // const us& Neq=(vvertex[0]).Neq;
    const us& Ns=gc->Ns;
    
    for(us k=0; k<vvertex.size();k++)
      {
	vvertex[k]->setRes(res.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1));
      }
  }
  
}		 // Namespace segment
