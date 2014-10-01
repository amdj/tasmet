#pragma once
#ifndef _GLOBALCONF_H_
#define _GLOBALCONF_H_




#include <memory>
#include "vtypes.h"
#include "material.h"
#include <assert.h>
#define MAXNF (30)
#define MINOMG (1e-3)
#define MAXOMG (1e5)


namespace tasystem{
  SPOILNAMESPACE
  // This class can be copied safely

  class TaSystem;
  
  class Globalconf{
    const TaSystem* thesys=NULL;
    d Mass=0;			/* Fluid mass in the system (should remain constant) */
    d omg;		// The "base" frequency in rad/s
    bool driven=true;
    us Nf__;			// Number of frequencies to solve for
    us Ns__;			// Corresponding number of time samples

  public:
    d T0,p0,rho0;			/* Reference temperature and pressure (used to initialize a lot of variables. */
    d c0;		// Typical cross-sectional area,
				// finite volume size, speed of sound,
				// deltax of volume
    d kappa;			// Artificial viscosity tuning factor,
				// typically between 0.25 and 0.75
    gases::Gas gas;
  public:
    Globalconf(){}
    Globalconf(us Nf,d freq,string Gas="air",d T0=293.15,d p0=101325.0,d kappa=1.0,bool driven=true);
    const us& Nf() const {return Nf__;}
    const us& Ns() const {return Ns__;}    
    static Globalconf airSTP(us Nf,d freq,d kappa=1.0);
    ~Globalconf(){TRACE(-5,"~Globalconf()");}
    const TaSystem* getSys()const { return thesys;}
    void setSys(TaSystem* sys) {thesys=sys;}
    bool isDriven() const {return driven;}
    void setDriven(bool d) { driven=d;}
    string Gastype;
    d getomg() const {return omg;}
    d getfreq() const {return omg/2/number_pi;}
    vd omgvec;    
    void setNf(us);
    void setp0(d p0){this->p0=p0; this->rho0=gas.rho(T0,p0);}
    void set(us Nf,d freq);	// Set data for new frequency and
				// number of samples
    void setomg(d omg);
    void setfreq(d freq);
    dmat iDFT; //inverse discrete Fourier transform matrix
    dmat fDFT; //forward discrete Fourier transform matrix
    dmat DDTfd;//Derivative in frequency domain
    dmat DDTtd;//Derivative in time domain
    dmat ddt; //Derivative matrix only nonzero frequency components
    dmat iddt; //Inverse of derivative matrix only nonzero frequency components

    void setMass(d mass){Mass=mass;}
    d getMass() const {return Mass;}
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
