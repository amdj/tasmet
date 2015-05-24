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

  const us maxNf=100;              // Maximum number of frequencies
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
  // const int rho=1;
  // const int m=2;
  // const int T=3;
  // const int p=4;
  // const int Ts=5;
  // Number of variables
  const int nvars_reserve=7;
  const int neqs_reserve=7;
  #endif

  
} // namespace constants

#endif

// Variables and their names
enum class Varnr{
  none,                         // None
    rho,                        // Density
    m,                          // Mass flow (rho*U)
    T,                          // Temperature
    p,                          // Pressure
    Ts,                         // Temperature of the solid
    mH,                         // Enthalpy flow (Watts)
    U,                          // Volume flow (m^3/s)
    u,                          // Velocity (U/Sf)
    mu,                         // Momentum flux
    Q,                          // Heat flow
    Qs                 // Solid heat Flow
    };

#ifndef SWIG
const char* toString(Varnr v);  // Convert Varnr to string. Defined in
                                // constants.cpp
#endif

namespace segment{
  enum Pos{left=0,right=1};
}
namespace tasystem{
  typedef segment::Pos Pos;
}
namespace tube{
  typedef segment::Pos Pos;

  #ifndef SWIG
  inline const char* posWord(Pos p){
    if(p==Pos::left)
      return "left";
    else
      return "right";
  }
  #endif
  
  #ifndef SWIG
  
  enum EqType{
    Con=0,			// Continuity
    Mom=1,			// Momentum
    Ene=2,			// Energy-like equation
    Ise=3,			// Isentropic
    Sta=4,			// State
    Sol=5,			// SolidEnergy
    Mu_is_m_u=6,      // momentumflow is massflow_squared div
                    // density*cs_area
    BcEq,
  };
  #endif  


}                // namespace tube

#endif // CONSTANTS_H
