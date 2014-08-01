/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "laminarduct.h"
#include "tubevertex.h"

// Tried to keep the method definition a bit in order in which a
  // tube is created, including all its components. First a tube is
  // created, which has a geometry and a global
  // configuration. Moreover, the tube has gridpoints, "LaminarDuctVertex"
  // instants. Of these, a tube has gp of them, stored in a vector. In
  // each gridpoint, variables live, which represent the current
  // solution. Moreover, we have equations in each gridpoint. More
  // precisely, in the final solution the continuity, momentum, energy
  // and a suitable equation of state should hold.
namespace tube {
  LaminarDuct::LaminarDuct(Geom geom):Tube(geom){
    // Fill vector of gridpoints with data:
    TRACE(13,"LaminarDuct constructor()...");
    type="LaminarDuct";
  }
  LaminarDuct::LaminarDuct(const LaminarDuct& other):LaminarDuct(other.geom){}
  LaminarDuct& LaminarDuct::operator=(const LaminarDuct& other){
    TRACE(13,"LaminarDuct copy assignment");
    Tube::operator=(other);
    // drag(geom);
    WARN("Do not use assignment operators for tubes");
    return *this;
  }  
  void LaminarDuct::init(const Globalconf& gc){
    assert(eq.size()==0);
    Tube::init(gc);
  }
  vector<const TubeEquation*> LaminarDuct::getEq() const {
    vector<const TubeEquation*> eq;
    eq.push_back(&c);
    eq.push_back(&m);
    eq.push_back(&e);
    eq.push_back(&s);
    eq.push_back(&se);
    return eq;
  }



  LaminarDuct::~LaminarDuct(){
    TRACE(15,"~LaminarDuct()");
    cleanup();
  }
  void LaminarDuct::cleanup(){}

  
} /* namespace tube */

