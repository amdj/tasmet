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
#include "tubebcvertex.h"
#include "drag.h"
#include "heat.h"


namespace tube{
  SPOILNAMESPACE
  using namespace segment;
  typedef vector<const TubeEquation*> EqVec;  

  class Tube:public Seg {

    // unique_ptrs can be checked if they have content by if(uniq_ptr).
    std::unique_ptr<TubeBcVertex> bcLeft; 
    std::unique_ptr<TubeBcVertex> bcRight;   
  public:
      
    Tube(const Geom& geom);
    Tube(const Tube& othertube); // Copy constructor copies everything!
    Tube& operator=(const Tube& othertube); // And again, we copy everything.
    ~Tube();
    void show(us showvertices=0) const;
    void addBc(const TubeBcVertex& vertex);
    virtual us getNDofs() const;
    virtual us getNEqs() const;    
    virtual EqVec getEqs() const=0;    	// Some derived class needs to define the equation array
    virtual const DragResistance& getDragResistance() const=0;
    virtual const HeatSource& getHeatSource() const=0;
    virtual string getType() const final {return string("Tube");}
    virtual void init(const Globalconf& gc);
    virtual void setDofNrs(us firstdofnr);
    virtual void setEqNrs(us firstdofnr);    
    virtual d getCurrentMass() const;	// Obtain current mass in system
    virtual vd dmtotdx() const; // Derivative of current mass in
				    // system to all dofs.
    
    vd getResAt(us varnr,us freqnr) const; // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    vd getErrorAt(us varnr,us freqnr) const; // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    friend class TubeVertex;    
  protected:
    void cleanup();
  private:
    void copyTube(const Tube&);
    TubeVertex* leftTubeVertex() const;
    TubeVertex* rightTubeVertex() const;    
  };				// Tube class

  
} /* namespace tube */

#endif /* TUBE_H_ */






