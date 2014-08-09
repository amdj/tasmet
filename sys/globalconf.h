#pragma once
#ifndef _GLOBALCONF_H_
#define _GLOBALCONF_H_




#include <memory>
#include <vtypes.h>
#include <material.h>
#include <assert.h>
#define MAXNF (30)
#define MINOMG (1e-3)
#define MAXOMG (1e5)


namespace tasystem{
  SPOILNAMESPACE
  // This class can be copied safely


  
  class Globalconf{
  public:
    Globalconf(){}
    Globalconf(us Nf,d freq,string Gas="air",d T0=293.15,d p0=101325.0,d Mass=0.0,d kappa=1.0);
    static Globalconf airSTP(us Nf,d freq,d Mass=0.0,d kappa=1.0);
    ~Globalconf(){TRACE(-5,"~Globalconf()");}
    // Globalconf(const Globalconf& o): Globalconf(o.Nf,o.freq,o.Gastype,o.T0,o.p0,o.Mach,o.S0,o.dx,o.Mass,o.kappa){}

    string Gastype;
    us Nf;			// Number of frequencies to solve for
    us Ns;			// Corresponding number of time samples

    d freq;
    d omg;		// The "base" frequency in rad/s
    d c0;		// Typical cross-sectional area,
				// finite volume size, speed of sound,
				// deltax of volume
    d kappa;			// Artificial viscosity tuning factor,
				// typically between 0.25 and 0.75

    
    d T0,p0,rho0;			/* Reference temperature and pressure (used to initialize a lot of variables. */
    d Mass;			/* Fluid mass in the system (should remain constant) */
    vd omgvec;
    gases::Gas gas;
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
    void show() const;
    //    void setgas(string g){ gas(g);}


  protected:
    void updateiDFT();
    void updatefDFT();
    void updateiomg();

    d oldomg; //Previous omega
  private:

  }; /* Class Globalconf */
}    // namespace tasystem
#endif /* _GLOBALCONF_H_ */
