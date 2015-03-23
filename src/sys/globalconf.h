// globalconf.h
//
// Author: J.A. de Jong 
//
// Description:
// Global configuration options
//////////////////////////////////////////////////////////////////////
#pragma once

#ifndef _GLOBALCONF_H_
#define _GLOBALCONF_H_

#include "vtypes.h"
#include "gas.h"
#include <assert.h>


void setLogger(int loglevel);

namespace tasystem{

  const us MAXNF=30;
  const d MINOMG=1e-3;
  const d MAXOMG=1e5;

  const d pSTP=101325;
  const d TSTP=293.15;

  #ifndef SWIG
  SPOILNAMESPACE
  #endif
  class Globalconf{
    d Mass=0;			/* Fluid mass in the system (should remain constant) */
    d omg;		// The "base" frequency in rad/s
    bool driven=true;
    us Nf_;			// Number of frequencies to solve for
    us Ns_;			// Corresponding number of time samples
    d T0_,p0_;			/* Reference temperature and pressure (used to initialize a lot of variables. */
    gases::Gas gas_;

  public:
    Globalconf(us Nf,d freq,const string& gasstring,d T0,d p0);
    static Globalconf airSTP(us Nf,d freq);
    static Globalconf heliumSTP(us Nf,d freq);
    
    const us& Nf() const {return Nf_;}
    const us& Ns() const {return Ns_;}    

    ~Globalconf(){TRACE(-5,"~Globalconf()");}
    d getomg() const {return omg;}
    d getfreq() const {return omg/2/number_pi;}
    d c0() const {return gas_.cm(T0_);}
    d rho0() const {return gas_.rho(T0_,p0_);}
    d deltanu0() const{ return sqrt(2*gas_.mu(T0_)/(rho0()*omg));}
    d T0() const {return T0_;}
    d p0() const {return p0_;}
    vd omgvec;    
    void setNf(us);
    #ifndef SWIG
    void set(us Nf,d freq);	// Set data for new frequency and
    // number of samples
    #endif
    void setomg(d omg);
    void setfreq(d freq);
    dmat iDFT; //inverse discrete Fourier transform matrix
    dmat fDFT; //forward discrete Fourier transform matrix
    dmat DDTfd;//Derivative in frequency domain
    dmat DDTtd;//Derivative in time domain
    dmat ddt; //Derivative matrix only nonzero frequency components
    dmat iddt; //Inverse of derivative matrix only nonzero frequency components

    void setGas(const string& mat){gas_=gases::Gas(mat);}
    const gases::Gas& gas() const {return gas_;}    

    void setMass(d mass){Mass=mass;}
    d getMass() const {return Mass;}

    void show() const;

  protected:

    void updateiDFT();
    void updatefDFT();
    void updateiomg();

    d oldomg; //Previous omega
  private:

  }; /* Class Globalconf */
}    // namespace tasystem

#endif /* _GLOBALCONF_H_ */
//////////////////////////////////////////////////////////////////////
