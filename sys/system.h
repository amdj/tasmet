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
#include <vtypes.h>
#include "globalconf.h"
#include "segconnection.h"
#include "segbase.h"
#define MAXSEGS 30

namespace tasystem{
  SPOILNAMESPACE

  class taSystem;
  using segment::SegBase;

  class taSystem{
  private:
    sdmat curjac;		// Current Jacobian matrix
    
    vector<std::unique_ptr<SegBase> > segs;		// The stl library does not have fixed-size vectors, so we fall back to basic C arrays since we do not want to use boost.
    vector<SegConnection> segConnections;
    arma::uvec::fixed<MAXSEGS> segfirstdof; // Vector containing the number of the first column corresponding to the first vertex of segment segfirstcol(i)
    arma::uvec::fixed<MAXSEGS> segndofs;  // Vector containe
    bool hasInit=false;
  public:

    Globalconf gc;    
  private:
    us getNDofs();	// Compute DOFS in system, set     
    
  public:
    taSystem(const Globalconf& g);
    ~taSystem();
    taSystem(const taSystem& o);
    taSystem& operator=(const taSystem& other);
    virtual taSystem* copy() const {return new taSystem(*this);}
    void show(bool showvertices=false);
    us getNSegs() const {return segs.size();}
    void connectSegs(us seg1,us seg2,SegCoupling);
    // System with a
    // vector of segments
    // ############################## ACCESS METHODS
    vd error();			// Total error vector
    vd getRes();			// Extract result vector
    void setRes(vd resvec);	// Set result vector
    void addSeg(const SegBase& s);	// Add a segment to the
					// system. It creates a copy
					// and ads it to segs by emplace_back.

    // ############################## ACCESS METHODS

    // void delseg(us n); // Not yet implementen. Delete a segment from the system (we have to determine how elaborated the API has to be.)
    void setGc(const Globalconf& gc); // Reset globalconf configuration
    dmat jac();		// Return Jacobian matrix    
    SegBase* operator[](us i) const;    
    SegBase* getSeg(us i) const; // Easier for cython wrapping
    void init();
  private:
    // A vector of boundary conditions is required
    void copytaSystem(const taSystem& other);
    void checkInit();
    void cleanup();
    // void setnodes(us segnr,us nL,us nR);
    // friend void copysegs(taSystem& to,const taSystem& from);
  };				// class System
  
} // namespace tasystem

#endif /* _SYSTEM_H_ */
