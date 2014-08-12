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
  Tube::Tube(const Geom& geom):Seg(geom){
    // Fill vector of gridpoints with data:
    TRACE(13,"Tube constructor()...");
  }
  void Tube::cleanup(){
    TRACE(13,"Tube::cleanup()");
    bcLeft.reset();
    bcRight.reset();
    vvertex.clear();
  }
  Tube::Tube(const Tube& other):Seg(other){
    copyTube(other);
  }
  Tube& Tube::operator=(const Tube& other){
    TRACE(13,"Tube copy assignment");
    cleanup();
    Seg::operator=(other);
    copyTube(other);
    return *this;
  }
  void Tube::show(bool showvertices) const {
    TRACE(18,"Tube::show()");
    Seg::show(showvertices);
    if(bcLeft)
      cout << "Left side contains internal boundary condition of type " << bcLeft->getType() << ".\n";
    if(bcRight)
      cout << "Right side contains internal boundary condition of type " << bcRight->getType() << ".\n";
    
  }
  void Tube::copyTube(const Tube& other){
    TRACE(13,"Tube::copyTube()");
    cleanup();
    if(other.bcLeft){
      TRACE(12,"Other has bc connected left side");
      bcLeft.reset(other.bcLeft->copy());
    }
    if(other.bcRight)
      bcRight.reset(other.bcLeft->copy());
  }
  
  void Tube::addBc(const TubeBcVertex& bc){
    TRACE(14,"Tube::addBc(bc)");
    if(bc.connectPos()==connectpos::left)
      {
	TRACE(12,"Bc connected left side");
	bcLeft.reset(bc.copy());
	if(bcLeft)
	  TRACE(12,"bcLeft is now "<< bool(bcLeft));
      }
    else if(bc.connectPos()==connectpos::right)
      {
	TRACE(12,"Bc connected right side");
	bcRight.reset(bc.copy());
      }
    else
      {      
	cout << "WARNING: bconnectbc(): Bc  not understood!\n";
      }
  }
  us Tube::getNDofs() const {
    TRACE(14,"Tube::getNDofs()");
    if(gc!=NULL)
      return vvertex.size()*gc->Ns*Neq;
    else
      return 0;
  }
  TubeVertex* Tube::leftTubeVertex() const{
    TRACE(13,"Tube::leftTubeVertex()");
    if(bcLeft){
      TRACE(13,"Tube::leftTubeVertex() returning a boundary vertex");
      return static_cast<TubeVertex*>(bcLeft->copy());
    }
    else{
      TRACE(13,"Tube::leftTubeVertex() returning an ordinary vertex");
      return new TubeVertex();
    }
  }
  TubeVertex* Tube::rightTubeVertex() const{
    if(bcRight)
      return static_cast<TubeVertex*>(bcRight->copy());
    else
      return new TubeVertex();
  }


  void Tube::init(const tasystem::Globalconf& g){
    TRACE(13,"Tube::Init()");
    Seg::init(g);
    if(vvertex.size()==0){
      TRACE(100,"Filling vertices. Current size:"<<vvertex.size());
      vvertex.emplace_back(leftTubeVertex());
      for(us i=1;i<geom.nCells-1;i++)
    	vvertex.emplace_back(new TubeVertex());
      vvertex.emplace_back(rightTubeVertex());

      us nVertex=vvertex.size();    
      assert(nVertex==geom.nCells);
      // And initialize again.
      for(us i=0;i<vvertex.size();i++){
	TubeVertex* cvertex=static_cast<TubeVertex*>(vvertex[i].get());
	TRACE(13,"Starting intialization of Vertex "<< i);
	if(i<nVertex-1) cvertex->setRight(*vvertex[i+1].get());
	if(i>0) cvertex->setLeft(*vvertex[i-1].get());
	cvertex->initTubeVertex(i,*this);
      }	// for
    }	// if
    else{
      WARN("Tube already initialized!");
    }
  } // Tube::init(gc)
  
    

  
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


