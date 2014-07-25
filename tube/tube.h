/*
 * tube.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#pragma once
#ifndef TUBE_H_
#define TUBE_H_
#include "seg.h"
#include "drag.h"
#include <vtypes.h>
#include <math_common.h>
// #include <Eigen/Sparse>


namespace tube{
  SPOILNAMESPACE

  using namespace segment;
  

  class Tube:public Seg {
  public:
    LaminarDragResistance drag;
    
    Tube(Geom geom);
    Tube(const Tube& othertube); // Copy constructor copies everything!
    Tube& operator=(const Tube& othertube); // And again, we copy everything.
    ~Tube();
    virtual void Init(const Globalconf& gc);
    virtual Vertex* makeVertex(us i,const Globalconf& gc);
    vd GetResAt(us varnr,us freqnr); // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    friend class TubeVertex;
    friend class Continuity;
    friend class Momentum;
  private:
    void cleanup();
  };				// Tube class

} /* namespace tube */

#endif /* TUBE_H_ */






