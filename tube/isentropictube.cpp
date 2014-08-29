/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "isentropictube.h"
#include "tubevertex.h"

// Tried to keep the method definition a bit in order in which a
  // tube is created, including all its components. First a tube is
  // created, which has a geometry and a global
  // configuration. Moreover, the tube has gridpoints, "IsentropicTubeVertex"
  // instants. Of these, a tube has gp of them, stored in a vector. In
  // each gridpoint, variables live, which represent the current
  // solution. Moreover, we have equations in each gridpoint. More
  // precisely, in the final solution the continuity, momentum, energy
  // and a suitable equation of state should hold.
namespace tube {
  IsentropicTube::IsentropicTube(Geom geom):Tube(geom){
    // Fill vector of gridpoints with data:
    TRACE(13,"IsentropicTube constructor()...");
  }
  IsentropicTube::IsentropicTube(const IsentropicTube& other):Tube(other){
    TRACE(13,"IsentropicTube copy cc");
  }
  IsentropicTube& IsentropicTube::operator=(const IsentropicTube& other){
    TRACE(13,"IsentropicTube copy assignment");
    Tube::operator=(other);
    // drag(geom);
    WARN("Do not use assignment operators for tubes");
    return *this;
  }  
  void IsentropicTube::init(const Globalconf& gc){
    Tube::init(gc);
    c.init(*this);
    m.init(*this);
    is.init(*this);
    s.init(*this);
    se.init(*this);
  }
  vector<const TubeEquation*> IsentropicTube::getEqs() const {
    vector<const TubeEquation*> eq;
    eq.push_back(&c);
    eq.push_back(&m);
    eq.push_back(&is);
    eq.push_back(&s);
    eq.push_back(&se);
    return eq;
  }

  IsentropicTube::~IsentropicTube(){
    TRACE(15,"~IsentropicTube()");
    cleanup();
  }
  void IsentropicTube::cleanup(){}

  
} /* namespace tube */
