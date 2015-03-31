// isentropictube.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
// Tried to keep the method definition a bit in order in which a
// tube is created, including all its components. First a tube is
// created, which has a geometry and a global
// configuration. Moreover, the tube has gridpoints, "IsentropicCell"
// instants. Of these, a tube has gp of them, stored in a vector. In
// each gridpoint, variables live, which represent the current
// solution. Moreover, we have equations in each gridpoint. More
// precisely, in the final solution the continuity, momentum, energy
// and a suitable equation of state should hold.
//////////////////////////////////////////////////////////////////////
#include "isentropictube.h"
#include "cell.h"
#include "globalconf.h"
#include "geom.h"

namespace tube {
  using tasystem::Globalconf;
  using tasystem::TaSystem;

  IsentropicTube::IsentropicTube(const Geom& geom):Tube(geom){
    // Fill vector of gridpoints with data:
    TRACE(13,"IsentropicTube constructor()...");
  }
  IsentropicTube::IsentropicTube(const IsentropicTube& other):Tube(other){
    TRACE(23,"IsentropicTube copy cc");
  }
  void IsentropicTube::init(const TaSystem& sys){
    Tube::init(sys);
    for(auto cell=cells.begin();cell!=cells.end();cell++){
      Cell& ccell=**cell;
      ccell.setIsentropic();
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

//////////////////////////////////////////////////////////////////////
