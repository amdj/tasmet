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
    TRACE(13,"Tube constructor()...");
    type="Tube";
  }
  Vertex* Tube::makeVertex(us i,const Globalconf& g){
    TRACE(13,"Tube::makeVertex("<< i <<",gc)");    
    return new TubeVertex();}
  void Tube::cleanup(){
    TRACE(13,"Tube::cleanup()");
    vvertex.clear();
    Nvertex=0;
    Ndofs=0;
  }
  Tube::Tube(const Tube& other):Tube(other.geom){}
  Tube& Tube::operator=(const Tube& other){
    TRACE(13,"Tube copy assignment");
    cleanup();
    geom=other.geom;
    // drag(geom);
    WARN("Do not use assignment operators for tubes");
    return *this;
  }
  void Tube::Init(const tasystem::Globalconf& g){
    TRACE(13,"Tube::Init()");
    Seg::Init(g);
    for (us i=0;i<Nvertex;i++){
      vvertex[i]->T.set(0,g.T0);
      vvertex[i]->rho.set(0,g.gas.rho(g.T0,g.p0));
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
  }
  

} /* namespace tube */

