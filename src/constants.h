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
  
} // namespace constants

#endif

// Variables and their names
// Unfortunately to let the code compile with Swig v 2.0, strongly
// typed enums are not supported. Therefore this is a normal
// enumerated type and not an enum class.
enum class Varnr {
  none,                         // None
    rho,                        // Density
    m,                          // Mass flow (rho*U)
    T,                          // Temperature
    p,                          // Pressure
    Ts,                         // Temperature of the solid
    Tw,                         // Temperature of the solid wall
    mH,                         // Enthalpy flow (Watts)
    U,                          // Volume flow (m^3/s)
    u,                          // Velocity (U/Sf)
    mu,                         // Momentum flux
    Q,                          // Heat flow
    Qs,                 // Solid heat Flow
    F,                 // A mechanical domain force [N]
    x,                // A mechanical displacement [m]
    Z,                 // A mechanical impedance [N/m]
    mEkin	       // Kinetic energy flow (Watts)
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
namespace mech{
  typedef segment::Pos Pos;
}
namespace duct{
  typedef segment::Pos Pos;

  #ifndef SWIG
  inline const char* posWord(Pos p){
    if(p==Pos::left)
      return "left";
    else
      return "right";
  }
  
  enum EqType{
    Con,			// Continuity
    Mom,			// Momentum
    Ene,			// Energy-like equation
    Ise,			// Isentropic
    Sta,			// State
    Sol,			// SolidEnergy
    SolTwEq,			// Solid wall temperature equation
    Mu_is_m_u,      // momentumflow is massflow_squared div
                    // density*cs_area
    BcEqP,
    BcEqu,
    BcEqStateBc,
  };
  #endif  


}                // namespace duct

#endif // CONSTANTS_H
