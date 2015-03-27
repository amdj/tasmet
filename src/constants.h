// constants.h
//
// Author: J.A. de Jong 
//
// Description:
// Definition of important constants
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef CONSTANTS_H
#define CONSTANTS_H 1

#ifndef SWIG

namespace constants {
  
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

  const d pSTP=101325;
  const d TSTP=293.15;
  
} // namespace constants

#endif

namespace tube{

  #ifndef SWIG
  const int RHONR=0;
  const int UNR=1;
  const int TNR=2;
  const int PNR=3;
  const int TSNR=4;
  // Number of variables
  const int NVARS=5;
  #endif

  enum varnr{rho=RHONR, U=UNR, T=TNR, p=PNR, Ts=TSNR,
             rhoL, rhoR, UL, UR, TL, TR, pL, pR, TsL, TsR
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
