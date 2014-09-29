/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "laminarduct.h"
#include "tubevertex.h"
#include "tube.h"
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

  LaminarDuct::LaminarDuct(const Geom& geom):Tube(geom),laminardrag(*this){
    // Fill vector of gridpoints with data:
    TRACE(13,"LaminarDuct constructor()...");
  }
  LaminarDuct::LaminarDuct(const LaminarDuct& o):Tube(o),
						 laminardrag(o.laminardrag)	\
						,heat(o.heat)
  {
    TRACE(13,"LaminarDuct copy constructor()...");
  }
  LaminarDuct& LaminarDuct::operator=(const LaminarDuct& o){
    TRACE(13,"LaminarDuct copy assignment");
    Tube::operator=(o);
    // drag(geom);
    cleanup();
    laminardrag=o.laminardrag;
    heat=o.heat;    
    return *this;
  }  
  void LaminarDuct::init(const Globalconf& gc){
    Tube::init(gc);
  }
  LaminarDuct::~LaminarDuct(){
    TRACE(15,"~LaminarDuct()");
    cleanup();
  }
  void LaminarDuct::cleanup(){ Tube::cleanup();}

  
} /* namespace tube */

