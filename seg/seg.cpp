#include "seg.h"
#include "bcvertex.h"
namespace segment{
  
  void coupleSegs(Seg& seg1,Seg& seg2,SegCoupling coupling){
    us seg1size=seg1.vvertex.size();
    us seg2size=seg2.vvertex.size();
    assert(seg1size>0);
    assert(seg2size>0);
    if (coupling==tailhead){
      // Seg1 is coupled with its tail to Seg2's head
      TRACE(3,"Coupling seg1 with its tail to the head of seg2");
      seg1.setRight(seg2);
      seg2.setLeft(seg1);
      seg1.vvertex[seg1size-1]->right=seg2.vvertex[0].get();
      seg2.vvertex[0]->left=seg1.vvertex[seg1size-1].get();

    }
    else if(coupling==headtail){
      // Seg2 is coupled with its tail to Seg1's head
      TRACE(3,"Coupling seg1 with its head to the tail of seg2");
      seg1.setLeft(seg2);
      seg2.setRight(seg1);
    }
    else if(coupling==tailtail){
      seg1.setRight(seg2);
      seg2.setRight(seg1);
    }
    else {
      // Coupling is headhead
      seg1.setLeft(seg2);
      seg2.setRight(seg1);
    }
  } // coupleSegs()
  
  Seg::Seg(Geom geom):SegBase(geom){
    TRACE(13,"Seg::Seg(Geom)");
    Ndofs=0;
    type="Seg";
    for(us j=0;j<MAXCONNECT;j++){
      left.push_back(NULL);
      right.push_back(NULL);
      // left[j]=NULL; right=[j]=NULL;
    }
    // The Jacobian matrix is larger than the number of dofs for the
    // connection terms other segments
    // us& Ns=gc.Ns;
  }
  Seg::Seg(const Seg& other): SegBase(other){
    TRACE(13,"Seg copy constructor");
    this->gc=other.gc;
    this->Ndofs=other.Ndofs;
    this->left=other.left;
    this->right=other.right;
    this->nleft=other.nleft;
    this->nright=other.nright;
    // The number is not copied, though
  }
  void Seg::Init(const tasystem::Globalconf& gc){
    TRACE(13,"Seg::Init()");
    this->gc=&gc;
    Nvertex=geom.Ncells;
    assert(Nvertex>0);
    Ndofs=geom.Ncells*gc.Ns*Neq;
    for(us i=0;i<Nvertex;i++){
      if(i<Nvertex-1)
	vvertex[i]->right=vvertex[i+1].get();
      if(i>0)
	vvertex[i-1]->left=vvertex[i].get();
      TRACE(13,"Starting intialization of Vertex "<< i);
      vvertex[i]->Init(i,gc,geom);
    }      
  }
  void Seg::show(){
    geom.show();
    showVertices();
  }
  void Seg::showVertices(){
    for(us i=0;i<Nvertex;i++)
      vvertex[i]->show();
  }
  void Seg::setLeft(const Seg& Left){
    TRACE(13,"Seg::SetLeft()");
    left[nleft]=&Left;
    nleft++;
  }
  void Seg::setRight(const Seg& Right){
    TRACE(13,"Seg::SetRight()");
    right[nright]=&Right;
    nright++;
  }
  void Seg::setLeftbc(Vertex* v){ // The segment owns the bc from then on!
    TRACE(13,"Seg::setLeftbc()-----EMPTY!");
    // delete vvertex[0];
    // vvertex[0]=
    TRACE(13,"setLeftBc(): v:"<<v);
    vvertex[0].reset(v);
    TRACE(13,"Segnumber of bc: " << ((BcVertex*) vvertex[0].get())->segNumber());
    // const us& Nvertex=geom.Nvertex;
    assert(Nvertex>0);
    vvertex[0]->right=vvertex[1].get();
    vvertex[1]->left=vvertex[0].get();
  }

  void Seg::setRightbc(Vertex* v){
    TRACE(13,"Seg::setRighbc()-----EMPTY!");
    assert(Nvertex>0);
    vvertex[Nvertex-1].reset(v);
    // delete vvertex[Nvertex-1];
    // const us& Nvertex=geom.Nvertex;
    // assert(Nvertex>0);
    vvertex[Nvertex-2]->right=v;
    vvertex[Nvertex-1]->left=vvertex[Nvertex-2].get();
  }    
  
  dmat Seg::Jac(){			// Return Jacobian matrix of error operator
    // sdmat Seg::Jac(){			// Return Jacobian matrix of error operator    
    TRACE(8," Seg::Jac().. ");
    
    // TRACE(-1,"Nvertex:"<<Nvertex);
    // sdmat Jac(Nvertex*Neq*Ns,Nvertex*Neq*Ns);
    const us& Ns=gc->Ns;
    dmat Jacobian (Nvertex*Neq*Ns,(Nvertex+2)*Neq*Ns,fillwith::zeros);    
    if(Jacobian.size()==0)
      Jacobian=dmat(Nvertex*Neq*Ns,(Nvertex+2)*Neq*Ns);    
    dmat vJac(Neq*Ns,3*Neq*Ns,fillwith::zeros);
    us firstrow,firstcol,lastrow,lastcol;
    TRACE(8,"Filling Segment Jacobian matrix...");
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
	  firstcol=Neq*Ns;
	  lastcol=4*Neq*Ns-1;
	} else{
	  firstcol=0;	  // Fill in
	  lastcol=3*Neq*Ns-1;	  
	}
      }	// j==0
      else {			// Last vertex
	if(vvertex[j]->right==NULL){
	  firstcol=(Nvertex-2)*Neq*Ns;
	  lastcol=(Nvertex+1)*Neq*Ns-1;	  
	}
	else{
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
