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
#include "tube.h"
#include "bcvertex.h"
#include "globalconf.h"
#include "systemhelpers.h"
#include "segconnection.h"

#define MAXSEGS 30

namespace tasystem{
  SPOILNAMESPACE

  using segment::Seg;
  using segment::BcVertex;
  using tube::Tube;


  class TAsystem;
  using segment::Seg;

  class TAsystem{
  private:
    vector<std::unique_ptr<Seg> > segs;		// The stl library does not have fixed-size vectors, so we fall back to basic C arrays since we do not want to use boost.
    vector<std::unique_ptr<BcVertex> > bcvertices;
    Globalconf gc;
    vector<SegConnection> segConnections;
    arma::uvec::fixed<MAXSEGS> segfirstdof; // Vector containing the number of the first column corresponding to the first vertex of segment segfirstcol(i)
    arma::uvec::fixed<MAXSEGS> segndofs;  // Vector containe
    us Ndofs=0;
    bool hasInit=false;
    
  private:
    void computeNdofs();	// Compute DOFS in system, set     
    
    
  public:
    TAsystem(const Globalconf& g);
    ~TAsystem();
    TAsystem(const TAsystem& o);
    TAsystem& operator=(const TAsystem& other);
    void show(bool showvertices=false);
    us getNSegs() const {return segs.size();}
    us getNBc() const {return bcvertices.size();}
    void connectSegs(us seg1,us seg2,SegCoupling);
    // System with a
    // vector of segments
    // ############################## ACCESS METHODS
    vd error();			// Total error vector
    vd getRes();			// Extract result vector
    void setRes(vd resvec);	// Set result vector
    void addSeg(const Seg& s);	// Add a segment to the system. It creates a copy.
    void addBc(const BcVertex& vertex);
    // ############################## ACCESS METHODS
    
    BcVertex* getBc(us nr) const;
    // void delseg(us n); // Not yet implementen. Delete a segment from the system (we have to determine how elaborated the API has to be.)
    void setGc(const Globalconf& gc); // Reset globalconf configuration
    dmat jac();		// Return Jacobian matrix    
    Seg* operator[](us i) const;    
    Seg* getSeg(us i) { return (*this)[i];} // Easier for cython wrapping
    void init();
  private:
    // A vector of boundary conditions is required

    void checkInit();
    void cleanup();
    // void setnodes(us segnr,us nL,us nR);
    // friend void copysegs(TAsystem& to,const TAsystem& from);
  };				// class System
  
} // namespace tasystem

#endif /* _SYSTEM_H_ */
