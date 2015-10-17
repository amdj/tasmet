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
#include "constants.h"
#include <exception>


namespace tasystem{

  #ifndef SWIG
  SPOILNAMESPACE
  #endif
  #ifdef SWIG
  %catches(std::exception,...) Globalconf::Globalconf(us Nf,d freq,const string& gasstring,d T0,d p0);
  %catches(std::exception,...) Globalconf::airSTP(us Nf,d freq);
  %catches(std::exception,...) Globalconf::heliumSTP(us Nf,d freq);
  #endif


  class Globalconf{
    d omg;		// The "base" frequency in rad/s
    us Nf_;			// Number of frequencies to solve for
    d T0_,p0_;			/* Reference temperature and pressure (used to initialize a lot of variables. */
    dmat fDFT_,iDFT_,DDTfd_;
    gases::Gas gas_;
    // d Wfo_=0;			// First order 'upwind' factor. If
				// Wfo=-1, interpolation is done from
				// the left side. If Wfo=0,
				// interpolation is second order. If
				// Wfo=1, interpolation is done from
				// the right side
  public:
    Globalconf(us Nf,d freq,const string& gasstring,d T0,d p0);
    static Globalconf airSTP(us Nf,d freq);
    static Globalconf heliumSTP(us Nf,d freq);
    
    us Nf() const {return Nf_;}
    us Ns() const {return 2*Nf_+1;}    

    ~Globalconf(){TRACE(-5,"~Globalconf()");}
    d getomg() const {return omg;}
    d getfreq() const {return omg/2/number_pi;}
    d meshPeclet(d dx,d u) const {return u*dx*rho0()*gas().cp(T0())/gas().kappa(T0());}
    d c0() const {return gas_.cm(T0_);}
    d rho0() const {return gas_.rho(T0_,p0_);}
    d deltanu0() const{ return sqrt(2*gas_.mu(T0_)/(rho0()*omg));}
    d deltanu0min() const{ return deltanu0()/sqrt((d) (max<us>(1,Nf())));}
    d T0() const {return T0_;}
    d p0() const {return p0_;}

    #ifndef SWIG
    void set(us Nf,d freq);	// Set data for new frequency and
    // number of samples
    #endif
    void setNf(us Nf);
    void setomg(d omg);
    void setfreq(d freq);
    const dmat& iDFT() const {return iDFT_;} //inverse discrete Fourier transform matrix
    const dmat& fDFT() const {return fDFT_;} //forward discrete Fourier transform matrix
    const dmat& DDTfd() const {return DDTfd_;}//Derivative in frequency domain

    void setGas(const string& mat){gas_=gases::Gas(mat);}
    const gases::Gas& gas() const {return gas_;}    

    void show() const;

  protected:
    void update();
    d oldomg; //Previous omega
  private:
  }; /* Class Globalconf */
}    // namespace tasystem

#endif /* _GLOBALCONF_H_ */
//////////////////////////////////////////////////////////////////////
