#pragma once
#include <vtypes.h>
#include <material.h>

namespace tasystem{

  class Globalconf{
  public:
    Globalconf(us Nf,d freq);
    Globalconf(us Nf,d freq,string Gas);
    Globalconf(us Nf,d freq,string Gas,d T0,d p0,d Mass);
    ~Globalconf();
    us Nf;			// Number of frequencies
    us Ns;			// Corresponding number of time samples
    gases::Gas gas;
    double freq;
    double omg;		// The "base" frequency in rad/s

    d T0,p0;			/* Reference temperature and pressure (used to initialize a lot of variables. */
    d Mass;			/* Fluid mass in the system (should remain constant) */


    void setMass(d mass){ Mass=mass;}
    void setfreq(d freq){freq=freq; omg=2*pi*freq; }
    void setp0(d p) { p0=p;}
    void setgas(gases::Gas g){ gas=g;}
    void show(){ }
    //    void setgas(string g){ gas(g);}

  }; /* Class Globalconf */
}    // namespace tasystem
