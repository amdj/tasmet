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
  
  using namespace segment;


  
  class Tube:public Seg {
  public:
    Tube(tasystem::Globalconf& g,Geom geom);
    Tube(const Tube& othertube); // Copy constructor to copy vertex
				 // vector ofpointers
    ~Tube();

    void setLeftbc(vertexptr v); // Set left boundary condition vertex
    void setRightbc(vertexptr v); // Set left boundary condition vertex    
    void setLeftbc(Vertex* v); // Set left boundary condition vertex, takes over ownership of object
    void setRightbc(Vertex* v); // Set right boundary condition vertex, takes over ownership of object
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






