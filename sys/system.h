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


#include <vtypes.h>
#include "tube.h"
#include "globalconf.h"


namespace tasystem{
  using segment::Seg;
  using arma::sp_mat;
  using math_common::esdmat;
  
  class TAsystem{
  public:
    TAsystem(Globalconf& g); // Initialize a
						 // System with a
						 // vector of segments
    vd Error();			// Total error vector
    vd GetRes();			// Extract result vector
    void SetRes(vd resvec);	// Set result vector
    void addseg(Seg& s);	// Add a segment to the system
    void delseg(us n); // Not yet implementen. Delete a segment from the system (we have to determine how elaborated the API has to be.)

    dmat Jac();		// Return Jacobian matrix    
    Seg& operator[](us i);    

  private:
    // A vector of boundary conditions is required
    void setnodes(us segnr,us nL,us nR);
    vector<segment::Seg*> segs;
    const Globalconf& gc;
    const us& Ns;

    // vector<us> startdof;	// Vector containing the starting degree of freedom for segment number # 
    // vector<us> enddof;		// Vector containing the last dof belonging to segment number #

    us Nsegs;			// Number of segments
    us Ndofs;


  };				// class System
  
} // namespace tasystem

#endif /* _SYSTEM_H_ */
