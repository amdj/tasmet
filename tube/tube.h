/*
 * tube.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#pragma once
#ifndef TUBE_H_
#define TUBE_H_
#include "globalconf.h"
#include "seg.h"
#include "geom.h"
#include "drag.h"
#include "tubevertex.h"
#include "var.h"
#include <vtypes.h>
#include <material.h>
#include <math_common.h>
// #include <Eigen/Sparse>




namespace tube{
  SPOILNAMESPACE
  using arma::sp_mat;  
  using namespace segment;

  class Tube:public Seg {
  public:
    Tube(Geom geom);
    Tube(const Tube& othertube); // Copy constructor only copies geometry.
    ~Tube();
    vd GetResAt(us varnr,us freqnr); // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
  protected:
    void Init(const tasystem::Globalconf& gc);
    LaminarDragResistance drag;
    friend class TubeVertex;
    friend class Continuity;
    friend class Momentum;
  private:

  };				// Tube class

} /* namespace tube */

#endif /* TUBE_H_ */






