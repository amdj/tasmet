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
#include "tubeequation.h"
#include "tubevertex.h"
#include "drag.h"
#include "heat.h"

namespace tube{
  SPOILNAMESPACE

  using namespace segment;
  typedef vector<const TubeEquation*> EqVec;  

  class Tube:public Seg {
  public:

      
    Tube(const Geom& geom);
    Tube(const Tube& othertube); // Copy constructor copies everything!
    Tube& operator=(const Tube& othertube); // And again, we copy everything.
    ~Tube();
    virtual EqVec getEq() const=0;    	// Some derived class needs to define the equation array
    virtual const DragResistance& getDragResistance() const=0;
    virtual const HeatSource& getHeatSource() const=0;
    virtual void setLeftBc(Vertex* v);
    virtual void setRightBc(Vertex* v);    
    virtual void init(const Globalconf& gc);
    vd getResAt(us varnr,us freqnr) const; // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
  protected:
    void cleanup();

  };				// Tube class

  
} /* namespace tube */

#endif /* TUBE_H_ */






