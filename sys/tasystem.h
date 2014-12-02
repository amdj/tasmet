// system.h, created March 19th 2014
// Author: J.A. de Jong

// A system contains one or more (un)connected (thermoacoustic)
// segments. Some segments require boundary conditions. A System()
// instance combines segments and boundary conditions, such that a
// full error vector can be created, a right hand side and the
// Jacobian.

#pragma once
#ifndef _SYSTEM_H_
#define _SYSTEM_H_
#define MAXNDOFS 600000
#include <memory>
#include "vtypes.h"
#include <boost/archive/text_oarchive.hpp>
#include "globalconf.h"
#include "segconnection.h"
#define MAXSEGS 30
#include "arma_eigen.h"

namespace tasystem{
  SPOILNAMESPACE

  typedef Eigen::VectorXd evd;
  typedef Eigen::SparseMatrix<double> esdmat;
  
  class TripletList;
  using segment::Seg;

  class TaSystem{
  private:
    vector<SegConnection> segConnections;
    bool hasInit=false;
  protected:
    vector<Seg*> segs;		
  public:
    Globalconf gc;    
  public:
    TaSystem(){}
    TaSystem(const Globalconf& g);
    TaSystem(const TaSystem& o);
    TaSystem& operator=(const TaSystem& other);
    virtual ~TaSystem();
    template <typename Archive>
    void serialize(Archive& ar,const unsigned int version){
      ar & gc;
    }
    virtual TaSystem* copy() const {return new TaSystem(*this);}
    virtual void show(us showvertices);
    virtual evd error();			// Total error vector
    virtual evd getRes();			// Extract result vector
    virtual void setRes(const vd& resvec);	// Set result vector
    virtual esdmat jac(d dummy=-1);		// Return Jacobian matrix
    virtual void init();
    us getNSegs() const {return segs.size();}
    void connectSegs(us seg1,us seg2,SegCoupling);
    void showJac(bool force=false);
    // System with a
    // vector of segments
    // ############################## ACCESS METHODS
    void setRes(const evd& res);
    void setNf(us);
    void addSeg(const Seg& s);	// Add a segment to the
					// system. It creates a copy
					// and ads it to segs by emplace_back.
    void addSeg(const std::vector<Seg*>&);
    void resetHarmonics();
    // ############################## ACCESS METHODS

    // void delseg(us n); // Not yet implementen. Delete a segment from the system (we have to determine how elaborated the API has to be.)
    void setGc(const Globalconf& gc); // Reset globalconf configuration

    Seg* operator[](us i) const;    
    Seg* getSeg(us i) const; // Easier for cython wrapping

    void setRes(const TaSystem& o);
    
    d getCurrentMass();	// Return current mass in system [kg]
    void checkInit(){		// Often called simple method: inline
      if(!hasInit){ init(); hasInit=true; }
    }
  protected:
    us getNDofs() const;	// Compute DOFS in system, set     
    us getNEqs() const;    
    void jacTriplets(TripletList&);
    void setDofEqNrs();
    // A vector of boundary conditions is required
    void copyTaSystem(const TaSystem& other);
    void cleanup();

  };				// class System
  
} // namespace tasystem

#endif /* _SYSTEM_H_ */
