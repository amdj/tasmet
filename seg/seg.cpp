#include "seg.h"
#include "bcvertex.h"
namespace segment{
  
  
  Seg::Seg(Geom geom):SegBase(geom){
    TRACE(13,"Seg::Seg(Geom)");
    Ndofs=0;
    type="Seg";
    // The Jacobian matrix is larger than the number of dofs for the
    // connection terms other segments
    // us& Ns=gc.Ns;
  }
  Seg::Seg(const Seg& other): SegBase(other){
    TRACE(13,"Seg copy constructor");
    this->gc=other.gc;
  }
  void Seg::Init(const tasystem::Globalconf& gc){
    TRACE(13,"Seg::Init()");
    // NO Do not clear segments! Boundary conditions have been added
    // vvertex.clear();
    Nvertex=geom.Ncells;
    Ndofs=geom.Ncells*gc.Ns*Neq;
    this->gc=&gc;

    if(vvertex.size()==0){
      for(us i=0;i<Nvertex;i++)
	vvertex.emplace_back(makeVertex(i,gc));
    }
    // And initialize again.
    for(us i=0;i<Nvertex;i++){
      TRACE(13,"Starting intialization of Vertex "<< i);
      if(i<Nvertex-1) vvertex[i]->right=vvertex[i+1].get();
      if(i>0) vvertex[i]->left=vvertex[i-1].get();
      vvertex[i]->Init(i,*this);
    }
  } // Seg::Init
  Vertex* Seg::makeVertex(us i,const Globalconf& gc){
    TRACE(13,"Seg::makeVertex()");
    return new Vertex();}
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
    for(us i=0;i<Nvertex;i++)
      vvertex[i]->show();
  }

  void Seg::setLeftbc(Vertex* v){ // The segment owns the bc from then on!
    TRACE(13,"Seg::setLeftbc()-----EMPTY!");
    assert(Nvertex>0);
    vvertex[0].reset(v);
    vvertex[0]->right=vvertex[1].get();
    vvertex[1]->left=vvertex[0].get();
    vvertex[0]->left=NULL;
    vvertex[0]->Init(0,*this);
  }

  void Seg::setRightbc(Vertex* v){
    TRACE(13,"Seg::setRighbc()-----EMPTY!");
    assert(Nvertex>0);
    vvertex[Nvertex-1].reset(v);
    vvertex[Nvertex-2]->right=v;
    vvertex[Nvertex-1]->left=vvertex[Nvertex-2].get();
    vvertex[Nvertex-1]->right=NULL;
    vvertex[Nvertex-1]->Init(Nvertex-1,*this);
  }    
  
  dmat Seg::Jac(){			// Return Jacobian matrix of error operator
    // sdmat Seg::Jac(){			// Return Jacobian matrix of error operator    
    TRACE(8," Seg::Jac() for Segment "<< getNumber() << ".");
    // TRACE(-1,"Nvertex:"<<Nvertex);
    // sdmat Jac(Nvertex*Neq*Ns,Nvertex*Neq*Ns);
    const us& Ns=gc->Ns;

    dmat Jacobian (Nvertex*Neq*Ns,(Nvertex+2)*Neq*Ns,fillwith::zeros);    
    if(Jacobian.size()==0)
      Jacobian=dmat(Nvertex*Neq*Ns,(Nvertex+2)*Neq*Ns);    
    dmat vJac(Neq*Ns,3*Neq*Ns,fillwith::zeros);
    us firstrow,firstcol,lastrow,lastcol;
    TRACE(8,"Filling Segment Jacobian matrix for segment "<< getNumber() <<"...");
    // #pragma omp parallel for
    for(us j=0;j<Nvertex;j++){			   // Fill the Jacobian
      TRACE(3,"Obtaining vertex Jacobian...");
      vJac=vvertex[j]->Jac();
      // The row height of a vertex jacobian matrix is Neq*Ns, The
      // column with is 3*Neq*Ns, since the neigbouring vertex has to
      // be found
      firstrow=j*Neq*Ns;
      lastrow=(j+1)*Neq*Ns-1;
      TRACE(3,"j:"<<j);
      if(j>0 && j<Nvertex-1){
	TRACE(0,"interior vertex jacobian, j="<<j);
	firstcol=(j)*Neq*Ns;
	lastcol=(j+3)*Neq*Ns-1;
      }
      else if(j==0){
	if(vvertex[j]->left==NULL){	// First node
	  TRACE(100,"First node not connected to other segments")
	  firstcol=Neq*Ns;
	  lastcol=4*Neq*Ns-1;
	} else{
	  TRACE(100,"First node IS connected to other segment")
	  firstcol=0;	  // Fill in
	  lastcol=3*Neq*Ns-1;	  
	}
      }	// j==0
      else {			// Last vertex
	if(vvertex[j]->right==NULL){
	  TRACE(100,"Last node not connected to other segments")
	  firstcol=(Nvertex-2)*Neq*Ns;
	  lastcol=(Nvertex+1)*Neq*Ns-1;	  
	}
	else{
	  TRACE(100,"Last node IS connected to other segment")
	  firstcol=(Nvertex-1)*Neq*Ns;
	  lastcol=Jacobian.n_cols-1;
	}
      }	// j==Nvertex-1
      TRACE(3,"Filling segment Jacobian with vertex subpart...");
      // cout << firstrow << " ";
      // cout << firstcol << " ";
      // cout << lastrow << " ";
      // cout << lastcol << " \n";
      // cout << "Jacsize:" << Jacobian.n_rows << " " << Jacobian.n_cols <<"\n";
      if(j==Nvertex-1)
	cout << "last vertex jac:\n"<< vJac;
      if(j==0)
	cout << "First vertex jac:\n"<< vJac;

      Jacobian.submat(firstrow,firstcol,lastrow,lastcol)=vJac;      
    }	// end for
    TRACE(8,"Segment Jacobian done.");
    return Jacobian;
  }
  vd Seg::GetRes(){
    TRACE(8,"Seg::GetRes()");

    vd Result(Ndofs,fillwith::zeros);
    us Ns=gc->Ns;
    for(us k=0; k<Nvertex;k++)
      {
	Result.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1)=vvertex[k]->GetRes();
      }
    return Result;
  }

  vd Seg::Error(){
    TRACE(8,"Seg::Error()");
    const us& Ns=gc->Ns;
    vd error(Ndofs,fillwith::zeros);
    for(us k=0; k<Nvertex;k++)
      {
	error.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1)=vvertex[k]->Error();
      }
    return error;
  }
  void Seg::SetRes(vd res){
    TRACE(8,"Seg::SetRes()");
    // const us& Neq=(vvertex[0]).Neq;
    const us& Ns=gc->Ns;
    for(us k=0; k<Nvertex;k++)
      {
	vvertex[k]->SetRes(res.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1));
      }
  }
  
}		 // Namespace segment
