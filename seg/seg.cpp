#include "seg.h"

namespace segment{
  static us totalnumber=0;
  
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
      }
  
  Seg::Seg(tasystem::Globalconf& g):gc(g),Ns(gc.Ns){
    TRACE(0,"Seg::Seg()");
    number=totalnumber;
    totalnumber++;
    Ndofs=Ncells=0;
    left=NULL; right=NULL;

    // The Jacobian matrix is larger than the number of dofs for the
    // connection terms other segments
    // us& Ns=gc.Ns;

  } // Seg constructor
  bool Seg::operator==(const Seg& other){return (this->number==other.number);}
  void Seg::setLeft(const Seg& Left){
    TRACE(0,"Seg::SetLeft()");
    left=&Left;
  }
  void Seg::setRight(const Seg& Right){right=&Right;}
  dmat Seg::Jac(){			// Return Jacobian matrix of error operator
    // sdmat Seg::Jac(){			// Return Jacobian matrix of error operator    
    TRACE(0," Seg::Jac().. ");
    // TRACE(-1,"Ncells:"<<Ncells);
    // sdmat Jac(Ncells*Neq*Ns,Ncells*Neq*Ns);
    dmat Jacobian (Ncells*Neq*Ns,(Ncells+2)*Neq*Ns,fillwith::zeros);    
    if(Jacobian.size()==0)
      Jacobian=dmat(Ncells*Neq*Ns,(Ncells+2)*Neq*Ns);    
    dmat vJac(Neq*Ns,3*Neq*Ns,fillwith::zeros);
    us firstrow,firstcol,lastrow,lastcol;
    TRACE(9,"Filling Segment Jacobian matrix...");
    // #pragma omp parallel for
    for(us j=0;j<Ncells;j++){			   // Fill the Jacobian
      TRACE(8,"Obtaining vertex Jacobian...");
      vJac=vvertex[j]->Jac();
      // The row height of a vertex jacobian matrix is Neq*Ns, The
      // column with is 3*Neq*Ns, since the neigbouring vertex has to
      // be found
      firstrow=j*Neq*Ns;
      lastrow=(j+1)*Neq*Ns-1;
      TRACE(8,"j:"<<j);
      if(j>0 && j<Ncells-1){
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
	  firstcol=(Ncells-2)*Neq*Ns;
	  lastcol=(Ncells+1)*Neq*Ns-1;	  
	}
	else{
	  firstcol=(Ncells-1)*Neq*Ns;
	  lastcol=Jacobian.n_cols-1;
	}
      }	// j==Ncells-1
      TRACE(8,"Filling segment Jacobian with vertex subpart...");
      // cout << firstrow << " ";
      // cout << firstcol << " ";
      // cout << lastrow << " ";
      // cout << lastcol << " \n";
      // cout << "Jacsize:" << Jacobian.n_rows << " " << Jacobian.n_cols <<"\n";
      Jacobian.submat(firstrow,firstcol,lastrow,lastcol)=vJac;      
    }	// end for
    TRACE(9,"Segment Jacobian done.");
    return Jacobian;
  }
  vd Seg::GetRes(){
    TRACE(0,"Seg::Get()");
    // const us& Neq=TubeVertex::Neq;
    const us& Ns=gc.Ns;
    vd Result(Ndofs,fillwith::zeros);
    for(us k=0; k<Ncells;k++)
      {
	Result.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1)=vvertex[k]->GetRes();
      }
    return Result;
  }

  vd Seg::Error(){
    TRACE(0,"Seg::Error()");
    const us& Ns=gc.Ns;
    vd error(Ndofs,fillwith::zeros);
    for(us k=0; k<Ncells;k++)
      {
	error.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1)=vvertex[k]->Error();
      }
    return error;
  }
  void Seg::SetRes(vd res){
    TRACE(0,"Seg::SetRes()");
    // const us& Neq=(vvertex[0]).Neq;
    const us& Ns=gc.Ns;
    for(us k=0; k<Ncells;k++)
      {
	vvertex[k]->SetRes(res.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1));
      }
  }
  
}		 // Namespace segment
