/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "tube.h"

  // Tried to keep the method definition a bit in order in which a
  // tube is created, including all its components. First a tube is
  // created, which has a geometry and a global
  // configuration. Moreover, the tube has gridpoints, "TubeVertex"
  // instants. Of these, a tube has gp of them, stored in a vector. In
  // each gridpoint, variables live, which represent the current
  // solution. Moreover, we have equations in each gridpoint. More
  // precisely, in the final solution the continuity, momentum, energy
  // and a suitable equation of state should hold.
namespace segment{
  Seg::Seg(tasystem::Globalconf& g):gc(g),vop(g.Nf,g.freq){
    TRACE(0,"Seg::Seg()");
    nL=0; nR=0;} // Seg constructor
}
namespace tube {
  Tube::Tube(tasystem::Globalconf& g,Geom geom):Seg(g),geom(geom),gas(gc.gas),drag(*this){
    // Fill vector of gridpoints with data:
    TRACE(5,"Tube constructor started, filling gridpoints vector...");

    Ncells=geom.Ncells;
    TRACE(0,"Ncells:"<<Ncells);
    vvertex=new Vertex*[Ncells];
    for(us i=0; i<Ncells;i++){
      TRACE(-1,"Tube vvertex i:"<<i);
      vvertex[i]=new TubeVertex(*this,i);
      // Link the array
      if(i>0){
	TRACE(-1,"Tube vvertex i-1:"<<i-1);
	vvertex[i]->left=vvertex[i-1];
	TRACE(-1,"Add right pointer to this one: " << vvertex[i]);
	vvertex[i-1]->right=vvertex[i];
      }

    }
    TRACE(5,"Tube constructor done");
    // globalconf instance is put in reference variable gc in
    // inherited class Seg
  }
  Tube::Tube(const Tube& o):Tube(o.gc,o.geom){
    TRACE(0,"Tube copy constructor");
    for(us i=0; i<Ncells;i++){
      *this->vvertex[i]=*o.vvertex[i];
    }
    // for this, we need to add to tubevertex copy constructor.
		// rule of Three
		// copy construcor
		// destructor
		// copy assignment operator
      // TODO fill this
  }
  void Tube::setLeftbc(TubeVertex* v){
    TRACE(0,"Tube::setLeftbc()");
    delete vvertex[0];
    vvertex[0]=v;
    vvertex[0]->right=vvertex[1];
    vvertex[1]->left=vvertex[0];
  }
  dmat Tube::Jac(){			// Return Jacobian matrix of error operator
    TRACE(0,"Tube::Jac()");
    us firstrow,firstcol,lastrow,lastcol;
    us& Ns=gc.Ns;
    dmat Jac(Ncells*Neq*Ns,Ncells*Neq*Ns,fillwith::zeros); // The Jacobian
						   // matrix is larger
						   // than the number
						   // of dofs for the
						   // connection terms
						   // to other
						   // segments
    dmat vJac;
    for(us j=0;j<Ncells;j++){			   // Fill the Jacobian
      vJac=vvertex[j]->Jac();
      // The row height of a vertex jacobian matrix is Neq*Ns, The
      // column with is 3*Neq*Ns, since the neigbouring vertex has to
      // be found
      firstrow=j*Neq*Ns;
      lastrow=(j+1)*Neq*Ns-1;

      if(j==0){
	firstcol=0;
	lastcol=2*Neq*Ns-1;
	Jac.submat(firstrow,firstcol,lastrow,lastcol)=vJac.submat(0,Neq*Ns,Neq*Ns-1,3*Neq*Ns-1);

      }
      else if(j==Ncells-1){
	firstcol=(Ncells-2)*Neq*Ns;
	lastcol=Jac.n_cols-1;//(Ncells-1)*Neq*Ns-1;
	TRACE(-1,firstrow<< " "<< firstcol << " " << lastrow << " " << lastcol);

	Jac.submat(firstrow,firstcol,lastrow,lastcol)=vJac.submat(0,0,Neq*Ns-1,2*Neq*Ns-1);


      }
      else{
	TRACE(-2,"interior vertex jacobian, j="<<j);
	firstcol=3*(j-1)*Neq*Ns;
	lastcol=3*j*Neq*Ns-1;

	TRACE(-1,"Jac rows"<<Jac.n_rows);
	TRACE(-1,"Jac ncols"<<Jac.n_cols);     	

	TRACE(0,firstrow<< " "<< firstcol << " " << lastrow << " " << lastcol);
	TRACE(-1,"vJac nrows"<<vJac.n_rows);   
	TRACE(-1,"vJac ncols"<<vJac.n_cols);   
	Jac.submat(firstrow,firstcol,lastrow,lastcol)=vJac;
	TRACE(-1,"submat nrows:"<< Jac.submat(firstrow,firstcol,lastrow,lastcol).n_rows);
	TRACE(-1,"submat ncols:"<< Jac.submat(firstrow,firstcol,lastrow,lastcol).n_cols);
      }
    }
    return Jac;
  }
  vd Tube::GetRes(){
    TRACE(0,"Tube::Get()");
    // const us& Neq=TubeVertex::Neq;
    const us& Ns=vop.Ns;
    vd Result(geom.Ncells*vop.Ns*Neq);
    for(us k=0; k<geom.Ncells;k++)
      {
	Result.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1)=vvertex[k]->GetRes();
      }
    return Result;
  }
  vd Tube::Error(){
    TRACE(0,"Tube::Error()");
    const us& Ns=vop.Ns;
    vd error(geom.Ncells*vop.Ns*Neq,fillwith::zeros);
    for(us k=0; k<geom.Ncells;k++)
      {
	error.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1)=vvertex[k]->Error();
      }
    return error;
  }
  void Tube::Set(vd res){
    TRACE(0,"Tube::Set");
    // const us& Neq=(vvertex[0]).Neq;
    const us& Ns=vop.Ns;
    for(us k=0; k<geom.Ncells;k++)
      {
	vvertex[k]->Set(res.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1));
      }
  }
  void Tube::Init(d T0,d p0){
    TRACE(0,"Tube::Init()");
    Vertex** v=vvertex;
    for (us i=0;i<Ncells;i++){
      TRACE(-1,"i:"<<i);
      (*v)->p.set(p0,0);
      (*v)->T.set(T0,0);
      (*v)->rho.set(gas.rho(T0,p0),0);
      v++;
    }
  }
  Tube::~Tube(){
    TRACE(-5,"Tube destructor started");
    if(Ncells>0){
      Vertex* v=vvertex[0];
      for (us i=0;i<Ncells;i++){
	TRACE(-6,"Deleting vertex..");
	delete vvertex[i];
      }
	TRACE(-6,"Deleting vertex array..");
      delete vvertex;  
    }

  }
  

} /* namespace tube */

