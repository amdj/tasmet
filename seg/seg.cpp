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
  void Seg::show(bool showVertices) const {
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
    if(showVertices==true)
      this->showVertices();
  }
  void Seg::showVertices() const {
    for(us i=0;i<vvertex.size();i++)
      vvertex[i]->show();
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
      assert(vJac.n_cols==3*Neq*Ns);      
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
      else if(j==nVertex-1){	// Last vertex
	if(getRight().size()==0){
	  TRACE(10,"Last node not connected to other segments")
	  firstcol=(nVertex-2)*Neq*Ns;
	  lastcol=(nVertex+1)*Neq*Ns-1;	  
	}
	else{
	  TRACE(10,"Last node IS connected to other segment")
	  firstcol=(j)*Neq*Ns;
	  lastcol=(j+3)*Neq*Ns-1;
	}
      }	// j==nVertex-1
      else
	{
	  WARN("Something went terribly wrong.");
	  exit(1);
	}
      
      TRACE(3,"Filling segment Jacobian with vertex subpart...");
      Jacobian.submat(firstrow,firstcol,lastrow,lastcol)=vJac;      
    }	// end for
    // cout <<"Segment" << getNumber() <<" Jacobian done. Jac is:\n"<< Jacobian;
    // cout << "Number of colums in this jacobian" << Jacobian.n_cols<<"\n";
    TRACE(8,"Segment Jacobian done.");
    return Jacobian;
  }
  vd Seg::getRes() const {
    TRACE(8,"Seg::GetRes()");
    vd Result(getNDofs(),fillwith::zeros);
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
    vd error(getNDofs(),fillwith::zeros);
    for(us k=0; k<nVertex;k++)
      {
	error.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1)=vvertex[k]->error();
      }
    return error;
  }
  vd Seg::domg() const{
    TRACE(8,"Seg::Error()");
    const us& Ns=gc->Ns;
    us nVertex=vvertex.size();    
    vd domg(getNDofs(),fillwith::zeros);
    for(us k=0; k<nVertex;k++)
      {
	domg.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1)=vvertex[k]->domg();
      }
    return domg;
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
