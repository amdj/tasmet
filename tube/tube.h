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
#include "varnr.h"

#include "drag.h"
#include "heat.h"


namespace tube{
  SPOILNAMESPACE
  using namespace segment;
  class TubeBcVertex;
  class DragResistance;
  class HeatSource;
  
  class Tube:public Seg {

    // unique_ptrs can be checked if they have content by if(uniq_ptr).
    std::unique_ptr<TubeBcVertex> bcLeft; 
    std::unique_ptr<TubeBcVertex> bcRight;   
  public:
      
    Tube(const Geom& geom);
    Tube(const Tube& othertube); // Copy constructor copies everything!
    Tube& operator=(const Tube& othertube); // And again, we copy everything.
    virtual ~Tube();
    void show(us showvertices=0) const;
    void addBc(const TubeBcVertex& vertex);
    virtual us getNDofs() const;
    virtual us getNEqs() const;    
    virtual const DragResistance& getDragResistance() const=0;
    virtual const HeatSource& getHeatSource() const=0;
    virtual string getType() const final {return string("Tube");}
    virtual void init(const Globalconf& gc);
    virtual void setDofNrs(us firstdofnr);
    virtual void setEqNrs(us firstdofnr);    
    virtual d getCurrentMass() const;	// Obtain current mass in system
    virtual void dmtotdx(vd&) const; // Derivative of current mass in
				    // system to all dofs.
    virtual void resetHarmonics();             // Set all higher
                                               // harmonic amplitudes
                                               // to zero
    virtual void setRes(const SegBase& other); // To copy from a
    virtual void updateNf();    
    // *similar* segment
    vd getResAt(us,us freqnr) const; // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    vd getResAt(varnr,us freqnr) const; // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    vd getErrorAt(us eqnr,us freqnr) const; // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    friend class TubeVertex;
    vd Htot() const;
    vd mtot() const;
    const TubeVertex& getTubeVertex(us i) const;
    vd interpolateResMid(varnr v,d x) const; // Amplitude data vectors
    vd interpolateResStaggered(varnr v,d x) const; // Amplitude data!!
    virtual vd dragCoefVec(us) const;              // return drag coefficient
  protected:
    void cleanup();
  private:
    void copyTube(const Tube&);
    TubeVertex* leftTubeVertex() const;
    TubeVertex* rightTubeVertex() const;    
  };				// Tube class

  
  
  inline Tube& asTube(SegBase& s){return dynamic_cast<Tube&>(s);}
  inline const Tube& asTube_const(const SegBase& s){return dynamic_cast<const Tube&>(s);}  
} /* namespace tube */

#endif /* TUBE_H_ */






