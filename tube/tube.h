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
#include "seg.h"
#include "geom.h"
#include "drag.h"
#include "vertex.h"
#include "../var/var.h"
#include <vtypes.h>
#include <material.h>
#include <math_common.h>

#include <Eigen/Sparse>




namespace tube{
  SPOILNAMESPACE
  using arma::sp_mat;  
  using segment::Seg;
  using segment::Vertex;
  class Tube:public Seg {
  public:
    Tube(tasystem::Globalconf& g,Geom geom);
    Tube(const Tube& othertube); // Copy constructor to copy vertex
				 // vector ofpointers
    ~Tube();

    void setLeftbc(TubeVertex* v); // Set left boundary condition vertex
    void setRightbc(TubeVertex* v); // Set left boundary condition vertex    
    Geom geom;			// The geometry
    gases::Gas& gas;		// The gas in the system. Reference variable to gc.gas
    vd GetResAt(us varnr,us freqnr); // Extrect a result vector for given variable number (rho,U,T,p,Ts) and frequency number.

  protected:
    void Init();
    LaminarDragResistance drag;
    friend class TubeVertex;
    friend class Continuity;
    friend class Momentum;
  private:

  };				// Tube class

} /* namespace tube */

#endif /* TUBE_H_ */






