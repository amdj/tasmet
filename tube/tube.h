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


namespace tube{
  SPOILNAMESPACE

  using namespace segment;
  

  class Tube:public Seg {
  public:
    vector<TubeEquation*> eq;
    
    Tube(Geom geom);
    Tube(const Tube& othertube); // Copy constructor copies everything!
    Tube& operator=(const Tube& othertube); // And again, we copy everything.
    ~Tube();
    virtual void setLeftBc(Vertex* v);
    virtual void setRightBc(Vertex* v);    
    virtual vector<const TubeEquation*> getEq() const=0;
    virtual void init(const Globalconf& gc);
    vd getResAt(us varnr,us freqnr) const; // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.

  private:
    void cleanup();
  };				// Tube class

  
} /* namespace tube */

#endif /* TUBE_H_ */






