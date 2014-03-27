// system.h, created March 19th 2014
// Author: J.A. de Jong

// A system contains one or more (un)connected (thermoacoustic)
// segments. Some segments require boundary conditions. A System()
// instance combines segments and boundary conditions, such that a
// full error vector can be created, a right hand side and the
// Jacobian.
#include <vtypes.h>
#include "tube/tube.h"
#include "globalconf.h"



namespace tasystem{
  using segment::Seg;
  class TAsystem{
  public:
    TAsystem(Globalconf& g); // Initialize a
						 // System with a
						 // vector of segments
    vd Error();			// Total error vector
    vd GetRes();			// Extract result vector
    void SetRes(vd resvec);	// Set result vector

    void addseg(Seg& s);
    void delseg(us n);

    dmat Jac();		// Return Jacobian matrix    
    Seg& operator[](us i);    
    
    ~TAsystem();
  protected:
    // A vector of boundary conditions is required
    vector<segment::Seg*> segs;
    vector<us> startdof;	// Vector containing the starting degree of freedom for segment number # 
    vector<us> enddof;		// Vector containing the last dof belonging to segment number #

    Globalconf& gc;
    us Nsegs;			// Number of segments
    us Ndofs;
    us& Ns;

  };				// class System
  
} // namespace tasystem

