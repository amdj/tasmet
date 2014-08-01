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
  Tube::Tube(Geom geom):Seg(geom){
    // Fill vector of gridpoints with data:
    TRACE(13,"Tube constructor()...");
    type="Tube";
  }
  void Tube::cleanup(){
    TRACE(13,"Tube::cleanup()");
    Seg::cleanup();
    vvertex.clear();
  }
  Tube::Tube(const Tube& other):Tube(other.geom){}
  Tube& Tube::operator=(const Tube& other){
    TRACE(13,"Tube copy assignment");
    geom=other.geom;
    cleanup();
    WARN("Do not use assignment operators for tubes");
    return *this;
  }
  void Tube::setLeftBc(Vertex* v){
    TRACE(13,"Tube::setLeftBc()");
    Seg::setLeftBc(v);
    us nVertex=vvertex.size();
    static_cast<TubeVertex*>(vvertex[0].get())->right=static_cast<TubeVertex*>(vvertex[1].get());
    static_cast<TubeVertex*>(vvertex[1].get())->left=static_cast<TubeVertex*>(vvertex[0].get());
    static_cast<TubeVertex*>(vvertex[0].get())->left=NULL;
    vvertex[0]->init(0,*this);
  }
  void Tube::setRightBc(Vertex* v){
    TRACE(13,"Tube::setRightBc()");
    Seg::setRightBc(v);
    us nVertex=vvertex.size();
    static_cast<TubeVertex*>(vvertex[nVertex-2].get())->right=static_cast<TubeVertex*>(vvertex[nVertex-1].get());
    static_cast<TubeVertex*>(vvertex[nVertex-1].get())->left=static_cast<TubeVertex*>(vvertex[nVertex-2].get());
    static_cast<TubeVertex*>(vvertex[nVertex-1].get())->right=NULL;
    vvertex[nVertex-1]->init(nVertex-1,*this);
  }
  void Tube::init(const tasystem::Globalconf& g){
    TRACE(13,"Tube::Init()");
    Seg::init(g);
    if(vvertex.size()==0){
      for(us i=0;i<geom.nCells;i++)
	vvertex.emplace_back(new TubeVertex());
    }
    us nVertex=vvertex.size();    
    // And initialize again.
    for(us i=0;i<vvertex.size();i++){
      TubeVertex* cvertex=static_cast<TubeVertex*>(vvertex[i].get());
      TRACE(13,"Starting intialization of Vertex "<< i);
      if(i<nVertex-1) cvertex->right=static_cast<TubeVertex*>(vvertex[i+1].get());
      if(i>0) cvertex->left=static_cast<TubeVertex*>(vvertex[i-1].get());
      cvertex->initTubeVertex(i,*this);
    }

  }
  vd Tube::getResAt(us varnr,us freqnr) const{
    const us& nCells=geom.nCells;
    vd res(nCells);
    assert(varnr<Neq);
    for(us i=0;i<nCells;i++){
      TubeVertex* cvertex=static_cast<TubeVertex*>(vvertex[i].get());
      res(i)=cvertex->vars[varnr]->operator()(freqnr);
    }
    return res;
  }
  
  Tube::~Tube(){
    TRACE(15,"~Tube()");
    cleanup();
  }
  

  
} /* namespace tube */


