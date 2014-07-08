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
#include <memory>
#include <vtypes.h>
#include "tube.h"
#include "globalconf.h"
#define MAXSEGS 30

namespace tasystem{
  SPOILNAMESPACE

  class TAsystem;
  typedef std::shared_ptr<tasystem::TAsystem> TAsystemptr;
  using segment::Seg;
  using segment::Segmentptr;

  class TAsystem{
  public:
    TAsystem(const Globalconf& g); // 
    ~TAsystem() {TRACE(-5,"~TAsystem()");}
    TAsystem(const TAsystem& o);
    // System with a
    // vector of segments
    vd Error();			// Total error vector
    vd GetRes();			// Extract result vector
    void SetRes(vd resvec);	// Set result vector
    void addseg(Segptr s);	// Add a segment to the system. From then on, the system owns the Segment. 
    void delseg(us n); // Not yet implementen. Delete a segment from the system (we have to determine how elaborated the API has to be.)
    void setGc(const Globalconf& gc); // Reset globalconf configuration
    dmat Jac();		// Return Jacobian matrix    
    Seg& operator[](us i);    
    arma::uvec segfirstcol=zeros<arma::uvec>(MAXSEGS);
    arma::uvec segndofs=zeros<arma::uvec>(MAXSEGS);

  private:
    // A vector of boundary conditions is required
    void setnodes(us segnr,us nL,us nR);
    vector<Segptr> segs;
    Globalconf gc;

    // vector<us> startdof;	// Vector containing the starting degree of freedom for segment number # 
    // vector<us> enddof;		// Vector containing the last dof belonging to segment number #

    us Nsegs;			// Number of segments
    us Ndofs;


  };				// class System
  
} // namespace tasystem

#endif /* _SYSTEM_H_ */
