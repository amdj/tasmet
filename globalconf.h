#pragma once
#include <vtypes.h>
#include <material.h>

namespace tasystem{
  SPOILNAMESPACE
  class Globalconf{
  public:
    Globalconf(us Nf,d freq,string Gas="air",d T0=293.15,d p0=101325.0,d Mach=1.0,d S0=1.0,d dx=1.0,d Mass=0.0,d kappa=0.25);
    ~Globalconf();
    gases::Gas gas;
    us Nf;			// Number of frequencies
    us Ns;			// Corresponding number of time samples

    d freq;
    d omg;		// The "base" frequency in rad/s
    d S0,V0,c0,dx,M;		// Typical cross-sectional area,
				// finite volume size, speed of sound,
				// deltax of volume
    d kappa;			// Artificial viscosity tuning factor,
				// typically between 0.25 and 0.75

    
    d T0,p0;			/* Reference temperature and pressure (used to initialize a lot of variables. */
    d Mass;			/* Fluid mass in the system (should remain constant) */
    vd omgvec;
    void set(us Nf,d freq);	// Set data for new frequency and number of samples
    dmat iDFT; //inverse discrete Fourier transform matrix
    dmat fDFT; //forward discrete Fourier transform matrix
    dmat DDTfd;//Derivative in frequency domain
    dmat DDTtd;//Derivative in time domain
    dmat ddt; //Derivative matrix only nonzero frequency components
    dmat iddt; //Inverse of derivative matrix only nonzero frequency components
    void setfreq(d freq);
    void setMass(d mass){ Mass=mass;}
    void setp0(d p) { p0=p;}
    void setgas(gases::Gas g){ gas=g;}
    void show(){ }
    //    void setgas(string g){ gas(g);}
  protected:
    void updateiDFT();
    void updatefDFT();
    void updateiomg();

    d oldomg; //Previous omega
  private:

  }; /* Class Globalconf */
}    // namespace tasystem
