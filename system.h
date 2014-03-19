// system.h, created March 19th 2014
// Author: J.A. de Jong

// A system contains one or more (un)connected (thermoacoustic)
// segments. Some segments require boundary conditions. A System()
// instance combines segments and boundary conditions, such that a
// full error vector can be created, a right hand side and the
// Jacobian.
#include <vtypes.h>
#include "tube/tube.h"

namespace system{
  class System{
  public:
    System(std::vector<segment::segment*> segs); // Initialize a
						 // System with a
						 // vector of segments

    ~System();
  protected:
    // A vector of boundary conditions is required
    us Nsegs;			// Number of segments
  };
} // namespace system
