// laminarduct.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////



#include "laminarduct.h"
#include "cell.h"
#include "tube.h"
#include "globalconf.h"
#include "geom.h"
#include "var.h"
// Tried to keep the method definition a bit in order in which a
  // tube is created, including all its components. First a tube is
  // created, which has a geometry and a global
  // configuration. Moreover, the tube has gridpoints, "LaminarDuctCell"
  // instants. Of these, a tube has gp of them, stored in a vector. In
  // each gridpoint, variables live, which represent the current
  // solution. Moreover, we have equations in each gridpoint. More
  // precisely, in the final solution the continuity, momentum, energy
  // and a suitable equation of state should hold.
using tasystem::Globalconf;
using tasystem::TaSystem;
using variable::var;

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

  void LaminarDuct::init(const TaSystem& sys){
    Tube::init(sys);
  }
  LaminarDuct::~LaminarDuct(){
    TRACE(15,"~LaminarDuct()");
  }
  vd LaminarDuct::dragCoefVec(us freqnr) const{
    TRACE(15,"LaminarDuct::drag_vec()");
    vd dragcoef(getNCells());
    var drag_varcoef(*gc);
    for(us i=0;i<dragcoef.size();i++){
      const Cell& cell=getCell(i);
      drag_varcoef.set(laminardrag.ComplexResistancecoef(cell));
      dragcoef(i)=drag_varcoef(freqnr);
    }
    return dragcoef;
  }
  
} /* namespace tube */

//////////////////////////////////////////////////////////////////////
