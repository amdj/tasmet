/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "tube.h"
#include "tubevertex.h"

// Tried to keep the method definition a bit in order in which a
  // tube is created, including all its components. First a tube is
  // created, which has a geometry and a global
  // configuration. Moreover, the tube has gridpoints, "TubeVertex"
  // instants. Of these, a tube has gp of them, stored in a vector. In
  // each gridpoint, variables live, which represent the current
  // solution. Moreover, we have equations in each gridpoint. More
  // precisely, in the final solution the continuity, momentum, energy
  // and a suitable equation of state should hold.
namespace tube {
  Tube::Tube(Geom geom):Seg(geom),drag(*this){
    // Fill vector of gridpoints with data:
    TRACE(13,"Tube constructor started, filling gridpoints vector...");

    // globalconf instance is put in reference variable gc in
    // inherited class Seg
    build();
    TRACE(13,"Tube constructor done");
  }
  void Tube::build(){
    TRACE(13,"Tube::build()");
    type="Tube";
    // vvertex=new Vertex*[Ncells];
    Nvertex=geom.Ncells;
    for(us i=0; i<Nvertex;i++){
      TRACE(12,"Creating Tube vvertex "<<i << " ...");
      vvertex.push_back(vertexptr(new TubeVertex()));
      // Link the array
      if(i>0){
	TRACE(-1,"Tube vvertex i-1:"<<i-1);
	vvertex[i]->left=vvertex[i-1].get();
	TRACE(-1,"Add right pointer to this one: " << vvertex[i]);
	vvertex[i-1]->right=vvertex[i].get();
      }
    }

  }
  void Tube::cleanup(){
    TRACE(13,"Tube::cleanup()");
    for(us i=0;i<Nvertex;i++)
      vvertex[i].reset();
    Nvertex=0;
    Ndofs=0;
  }
  
  Tube::Tube(const Tube& other):Tube(other.geom){
    TRACE(13,"Tube copy constructor");
    this->gc=other.gc;
    this->Ndofs=other.Ndofs;
    this->left=other.left;
    this->right=other.right;
    this->nleft=other.nleft;
    this->nright=other.nright;
    TRACE(13,"Tube copy constructor done.");
    // First runs the Seg copy constructor. This copies the pointers to left and righ segment of this one. Then the Seg base constructor is called from the Seg Cc. After that the Tube constructor is called to create the vertices.
  }
  Tube& Tube::operator=(const Tube& other){
    TRACE(0,"Tube copy assignment");
    cleanup();
    newgeom(other.geom);
    // drag(geom);
    build();
    return *this;
  }
  void Tube::Init(const tasystem::Globalconf& g){
    TRACE(13,"Tube::Init()");
    Seg::Init(g);
    for (us i=0;i<Nvertex;i++){
      // TRACE(-1,"i:"<<i);
      vvertex[i]->T.set(g.T0,0);
      vvertex[i]->rho.set(g.gas.rho(g.T0,g.p0),0);
    }
  }


  vd Tube::GetResAt(us varnr,us freqnr){
    const us& Ncells=geom.Ncells;
    vd res(Ncells);
    assert(varnr<Neq);
    for(us i=0;i<Ncells;i++){
      res(i)=vvertex[i]->vars[varnr]->operator()(freqnr);
    }
    return res;
  }
  Tube::~Tube(){
    TRACE(15,"~Tube()");
    cleanup();
    // if(Ncells>0){
      // Vertex* v=vvertex[0];
      // for (us i=0;i<Ncells;i++){
	// TRACE(-6,"Deleting vertex..");
	// delete vvertex[i];
      // }
	// TRACE(-6,"Deleting vertex array..");
      // delete vvertex;  
    // }

  }
  

} /* namespace tube */

