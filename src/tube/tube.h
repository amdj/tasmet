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

namespace tasystem{
  class Jacobian;
  class TaSystem;
}
namespace tube{
  SPOILNAMESPACE
  class DragResistance;
  class HeatSource;
  class TubeVertex;
  class Geom;

  class TubeBcVertex;

  class Tube:public segment::Seg {
    void showVertices(us detailnr) const ;   
    Geom* geom_;			// The geometry    
    Tube& operator=(const Tube&); // no copies allowed

  protected:
    std::vector<TubeVertex*> vvertex;
  public:
    Tube(const Geom& geom);
    Tube(const Tube& other); // Copy constructor copies everything!
    us getNCells() const;
						   // to some segment
						   // on right side
    virtual void init(const tasystem::TaSystem&);
    virtual ~Tube();

    void setRes(const Seg& other); // To copy from a
    void setResVar(varnr,us freqnr,d value);
    void setResVar(varnr,const vd& value);
    void show(us showvertices=0) const;
    const Geom& geom() const;
    vd getResAt(us,us freqnr) const; // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    vd getResAt(varnr,us freqnr) const; // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    vd getErrorAt(us eqnr,us freqnr) const; // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    const TubeVertex& operator[](us i) const;
    const TubeBcVertex& leftVertex() const;
    const TubeBcVertex& rightVertex() const;
    us getNVertex() const {return vvertex.size();}    
    vd Htot() const;
    vd mtot() const;
    const TubeVertex& getTubeVertex(us i) const;
    vd interpolateResMid(varnr v,d x) const; // Amplitude data vectors
    vd interpolateResStaggered(varnr v,d x) const; // Amplitude data!!

    virtual us getNDofs() const;
    virtual us getNEqs() const;    
    virtual const DragResistance& getDragResistance() const=0;
    virtual const HeatSource& getHeatSource() const=0;
    virtual string getType() const final {return "Tube";}
    virtual void setDofNrs(us firstdofnr);
    virtual void setEqNrs(us firstdofnr);    
    virtual d getCurrentMass() const;	// Obtain current mass in
                                        // system
    virtual void dmtotdx(vd&) const; // Derivative of current mass in
				    // system to all dofs.
    virtual void resetHarmonics();             // Set all higher
                                               // harmonic amplitudes
                                               // to zero
    virtual void domg(vd& tofill) const;
    virtual vd error() const;
    virtual void jac(tasystem::Jacobian& tofill) const;
    virtual vd getRes() const;
    virtual d getRes(us dofnr) const;
    virtual void setRes(const vd& res);
    virtual void updateNf();    

    // *similar* segment
    virtual vd dragCoefVec(us) const;              // return drag
                                                   // coefficient
  private:
    void cleanup_vvertex();
  };				// Tube class

  
  inline Tube& asTube(segment::Seg& s){return dynamic_cast<Tube&>(s);}
  inline const Tube& asTube_const(const segment::Seg& s){return dynamic_cast<const Tube&>(s);}  
} /* namespace tube */

#endif /* TUBE_H_ */






