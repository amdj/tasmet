// constants.h
//
// Author: J.A. de Jong 
//
// Description:
// Definition of important constants
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifndef SWIG

template<typename T>
T max(T t1,T t2) { return t1>t2?t1:t2;}
template<typename T>
T min(T t1,T t2) { return t1>t2?t2:t1;}


namespace constants {
  typedef double d;
  typedef unsigned us;

  const int mingp=4;              // Minimum number of gridpoints
  const int maxgp=3000;           // Maximum number of gridpoints

  const us maxNf=30;              // Maximum number of frequencies
  const d minomg=1e-3;            // Minimal oscillation frequency
  const d maxomg=1e5;

  const int maxsegs=30;           // Maximum number of segments in a TaSystem
  const int maxndofs=600000;      // Maximum number of DOFS

  const d minp=1e0;
  const d maxp=1e7;
  const d minT=2;                 // Minimal temperature
  const d maxT=2000;              // Maximal temperature

  const d p0=101325;
  const d T0=293.15;

  #ifndef SWIG
  // These variable numbers are important, as they determine the
  // position of these variables in the array in cell.h
  const int rho=0;
  const int rhoU=1;
  const int T=2;
  const int p=3;
  const int Ts=4;
  // Number of variables
  const int nvars=5;
  const int neqs=5;
  #endif

  
} // namespace constants

#endif

namespace segment{
  enum Pos{left=0,right=1};
}

namespace tube{
  typedef segment::Pos Pos;

  enum Varnr{rho=constants::rho,
             rhoU=constants::rhoU,
             T=constants::T,
             p=constants::p,
             Ts=constants::Ts,
             U
  };

  enum physquant{massFlow,
                 momentumFlow,
                 energyFlow,
                 heatFlow,
                 solidHeatFlow,
                 rhoRT,
  };
}                // namespace tube

#endif // CONSTANTS_H
