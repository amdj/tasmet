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
#include "constants.h"
#include "drag.h"
#include "heat.h"

namespace tasystem{
  class Jacobian;
  class TaSystem;
}

namespace tube{
  #ifndef SWIG
  SPOILNAMESPACE


  class DragResistance;
  class HeatSource;
  class TubeVertex;
  class Geom;

  class TubeBcVertex;
  #endif

  class Tube:public segment::Seg {
    void showVertices(us detailnr) const ;   
    Geom* geom_=nullptr;			// The geometry    


  protected:
    Tube(const Tube& other);
    std::vector<TubeVertex*> vvertex;
  public:
    Tube(const Geom& geom) throw(std::exception);

    Tube& operator=(const Tube&)=delete; // no copies allowed
    virtual ~Tube();          // Define this class as abstract

    const Geom& geom() const;
    vd Htot() const throw(std::exception);
    void setResVar(varnr,us i,us freqnr,d value);
    void setResVar(varnr,us freqnr,const vd& value);
    vd getValue(varnr,us freqnr) const throw(std::exception); // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    vc getValueC(varnr,us freqnr) const throw(std::exception); // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    vd getErrorAt(us eqnr,us freqnr) const throw(std::exception); // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.

    // Methods not exposed to swig
    #ifndef SWIG
    us getNCells() const;
    virtual bool init(const tasystem::TaSystem&);
    void setRes(const segment::Seg& other); // To copy from a
    void show(us showvertices=0) const;


    us getNVertex() const {return vvertex.size();}    
    vd interpolateResMid(varnr v,d x) const; // Amplitude data vectors
    vd interpolateResStaggered(varnr v,d x) const; // Amplitude data!!

    virtual us getNDofs() const;
    virtual us getNEqs() const;    
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

    virtual vd getRes() const;
    virtual d getRes(us dofnr) const;
    virtual void setRes(const vd& res);
    virtual void updateNf();    

    // *similar* segment
    virtual vd dragCoefVec(us) const;              // return drag
                                                   // coefficient

    const TubeVertex& operator[](us i) const;
    const TubeBcVertex& leftVertex() const;
    const TubeBcVertex& rightVertex() const;
    const TubeVertex& getTubeVertex(us i) const;
    virtual void jac(tasystem::Jacobian& tofill) const;
    virtual const DragResistance& getDragResistance() const=0;
    virtual const HeatSource& getHeatSource() const=0;
    #endif
  private:
    void cleanup_vvertex();
  };				// Tube class

  #ifndef SWIG  
  inline Tube& asTube(segment::Seg& s){return dynamic_cast<Tube&>(s);}
  inline const Tube& asTube_const(const segment::Seg& s){return dynamic_cast<const Tube&>(s);}
  #endif
} /* namespace tube */

#endif /* TUBE_H_ */






