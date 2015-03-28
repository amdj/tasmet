/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "isentropictube.h"
#include "tubevertex.h"
#include "globalconf.h"
#include "geom.h"
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
  using tasystem::Globalconf;
  using tasystem::TaSystem;

  IsentropicTube::IsentropicTube(const Geom& geom):Tube(geom){
    // Fill vector of gridpoints with data:
    TRACE(13,"IsentropicTube constructor()...");
  }
  IsentropicTube::IsentropicTube(const IsentropicTube& other):Tube(other){
    TRACE(13,"IsentropicTube copy cc");
  }
  void IsentropicTube::init(const TaSystem& sys){
    Tube::init(sys);
    for(auto vertex=vvertex.begin();vertex!=vvertex.end();vertex++){
      TubeVertex& cvertex=**vertex;
      cvertex.setIsentropic();
    }
    setInit(true);
  }
  vd IsentropicTube::dragCoefVec(us i) const {
    return zeros<vd>(geom().nCells());
  }

  IsentropicTube::~IsentropicTube(){
    TRACE(15,"~IsentropicTube()");
  }

  
} /* namespace tube */

