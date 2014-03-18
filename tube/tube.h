/*
 * tube.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#pragma once
#ifndef TUBE_H_
#define TUBE_H_
#include "../globalconf.h"
#include "geom.h"
#include "drag.h"
#include "vertex.h"
#include "../var/var.h"
#include <vtypes.h>
#include <material.h>
#include <math_common.h>

namespace segment{
  class Seg{
  public:
    Seg(globalconf::Globalconf& gc); // nL,nR initiated as 0
    ~Seg(){}
    us nL,nR;
    globalconf::Globalconf& gc;	// Global configuration of the system
    void setnodes(us n1,us n2){ nL=n1; nR=n2;}
  };

} // Namespace segment

namespace tube{
  class Tube:public segment::Seg {
  public:
    Tube(globalconf::Globalconf& g,Geom geom);
    ~Tube();
    void Init(d T0,d p0);
    variable::varoperations vop;
    Geom geom;			// The geometry
    gases::Gas& gas;		// The gas in the system. Reference variable to globalconf.gas
    vector<TubeVertex> vvertex; // Vector of vertices
    dmat Jacobian();
    vd Get();
    vd Error();
    void Set(vd res);
  protected:
    LaminarDragResistance drag;
    friend class TubeVertex;
    friend class Continuity;
    friend class Momentum;
  };				// Tube class

} /* namespace tube */

#endif /* TUBE_H_ */
