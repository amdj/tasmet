#include "common/vtypes.h"
#include "common/material.h"

namespace globalconf{

  class Globalconf{
  public:
    Globalconf(us Nf,d freq);
    ~Globalconf();
    double omega;		// The "base" frequency in rad/s
    double freq;
    us Nf;			// Number of frequencies
    us Ns;			// Corresponding number of time samples
    d T0,p0;			/* Reference temperature and pressure (used to initialize a lot of variables. */
    d Mass;			/* Fluid mass in the system (should remain constant) */
    gases::Gas& gas;
  }; /* Class Globalconf */
}
