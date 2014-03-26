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
    Seg(tasystem::Globalconf& gc); // nL,nR initiated as 0
    virtual ~Seg(){}
    us nL,nR;
    tasystem::Globalconf& gc;	// Global configuration of the system
    variable::varoperations vop;
    virtual vd Error()=0;
    virtual vd GetRes()=0;
    virtual dmat Jac()=0;
    void setnodes(us n1,us n2){ nL=n1; nR=n2;}
    Vertex** vvertex; // Vector of vertices
    us Ncells;
  };

} // Namespace segment

namespace tube{
  using segment::Seg;
  using segment::Vertex;
  class Tube:public Seg {
  public:
    Tube(tasystem::Globalconf& g,Geom geom);
    Tube(const Tube& othertube); // Copy constructor to copy vertex
				 // vector ofpointers
    ~Tube();
    void Init(d T0,d p0);

    void setLeftbc(TubeVertex* v); // Set left boundary condition vertex
    Geom geom;			// The geometry
    gases::Gas& gas;		// The gas in the system. Reference variable to globalconf.gas

    dmat Jac();
    vd GetRes();
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






