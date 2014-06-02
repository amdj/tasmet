/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "tube.h"
#include <math_common.h>

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
  Tube::Tube(const tasystem::Globalconf& g,Geom geom):Seg(g),geom(geom),gas(g.gas),drag(*this){
    // Fill vector of gridpoints with data:
    TRACE(5,"Tube constructor started, filling gridpoints vector...");
    type="Tube";
    Ncells=geom.Ncells;
    Ndofs=Ncells*gc.Ns*Neq;
    TRACE(0,"Ncells:"<<Ncells);
    // vvertex=new Vertex*[Ncells];
    for(us i=0; i<Ncells;i++){
      TRACE(-1,"Tube vvertex i:"<<i);
      vvertex.push_back(vertexptr(new TubeVertex(*this,i)));
      // Link the array
      if(i>0){
	TRACE(-1,"Tube vvertex i-1:"<<i-1);
	vvertex[i]->left=vvertex[i-1].get();
	TRACE(-1,"Add right pointer to this one: " << vvertex[i]);
	vvertex[i-1]->right=vvertex[i].get();
      }

    }
    TRACE(5,"Tube constructor done");
    // globalconf instance is put in reference variable gc in
    // inherited class Seg
    Init();
  }
  void Tube::Init(){
    TRACE(0,"Tube::Init()");
    // Vertex** v=vvertex;
    for (us i=0;i<Ncells;i++){
      TRACE(-1,"i:"<<i);
      vvertex[i]->T.set(gc.T0,0);
      vvertex[i]->rho.set(gas.rho(gc.T0,gc.p0),0);
      // v++;
    }
  }

  Tube::Tube(const Tube& o):Tube(o.gc,o.geom){
    TRACE(0,"Tube copy constructor");
    for(us i=0; i<Ncells;i++){
      // this->vvertex[i].reset(new TubeVertex(*(o.vvertex[i])));
    }
    // for this, we need to add to tubevertex copy constructor.
		// rule of Three
		// copy construcor
		// destructor
		// copy assignment operator
      // TODO fill this
  }

  vd Tube::GetResAt(us varnr,us freqnr){
    vd res(Ncells);
    assert(varnr<Neq);
    for(us i=0;i<Ncells;i++){
      res(i)=vvertex[i]->vars[varnr]->operator()(freqnr);
    }
    return res;
  }
  Tube::~Tube(){
    TRACE(-5,"Tube destructor started");
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

