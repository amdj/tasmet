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
#define MAXNDOFS 20000
#include <memory>
#include "vtypes.h"

#include "globalconf.h"
#include "segconnection.h"
#include "segbase.h"
#define MAXSEGS 30
#include "arma_eigen.h"

namespace tasystem{
  SPOILNAMESPACE

  typedef Eigen::VectorXd evd;
  typedef Eigen::SparseMatrix<double> esdmat;
  
  class TaSystem;
  using segment::SegBase;

  class TaSystem{
  private:
    vector<SegConnection> segConnections;
    arma::uvec::fixed<MAXSEGS> segfirstdof; // Vector containing the number of the first column corresponding to the first vertex of segment segfirstcol(i)
    arma::uvec::fixed<MAXSEGS> segndofs;  // Vector containe
    bool hasInit=false;
  protected:
    vector<std::unique_ptr<SegBase> > segs;		
  public:
    Globalconf gc;    
  protected:
    us getNDofs() const;	// Compute DOFS in system, set     
    
  public:
    TaSystem(const Globalconf& g);
    TaSystem(const TaSystem& o);
    TaSystem& operator=(const TaSystem& other);
    virtual ~TaSystem();
    virtual TaSystem* copy() const {return new TaSystem(*this);}
    virtual void show(bool showvertices);
    us getNSegs() const {return segs.size();}
    void connectSegs(us seg1,us seg2,SegCoupling);
    // System with a
    // vector of segments
    // ############################## ACCESS METHODS
    virtual evd error();			// Total error vector
    virtual evd getRes();			// Extract result vector
    virtual void setRes(const vd& resvec);	// Set result vector
    virtual esdmat jac();		// Return Jacobian matrix    

    
    void addSeg(const SegBase& s);	// Add a segment to the
					// system. It creates a copy
					// and ads it to segs by emplace_back.

    // ############################## ACCESS METHODS

    // void delseg(us n); // Not yet implementen. Delete a segment from the system (we have to determine how elaborated the API has to be.)
    void setGc(const Globalconf& gc); // Reset globalconf configuration

    SegBase* operator[](us i) const;    
    SegBase* getSeg(us i) const; // Easier for cython wrapping
    virtual void init();

    
    d getCurrentMass() const;	// Return current mass in system [kg]
    void checkInit(){		// Often called simple method: inline
      if(!hasInit){ init(); hasInit=true; }
    }
  private:
    // A vector of boundary conditions is required
    void copyTaSystem(const TaSystem& other);
    void cleanup();

  };				// class System
  
} // namespace tasystem

#endif /* _SYSTEM_H_ */
