// file: tubebc.h, created March 20th, 2014
// Author: J.A. de Jong

// bcvertex.h: external boundary conditions for tubes. This file
// contains the implementation of typical external boundary conditions
// for tubes as a custom vertex. Examples are adiabatic walls, isothermal walls and an
// adiabatic open pressure boundary conditions.

#include "tube.h"
#include "vertex.h"
#include "continuityeq.h"
#include "momentumeq.h"
#include "energyeq.h"

namespace tube{
  class LeftPressure:public TubeVertex
  {
  public:
    LeftPressure(const Tube& t);
    LeftPressure(const Tube&t,variable::var& pres);
    LeftPressure(const Tube&t,variable::var& pres,variable::var& temp);

    ~LeftPressure();
    void Init();
    // virtual vd cbcsource();
    virtual vd msource() const;	// Prescribed left pressure
    virtual vd esource() const;	// Same
    variable::var pL;			// Pressure boundary condition
    variable::var TL;			// Temperature boundary conditions
  };

  class RightImpedanceMomentumEq:public Momentum{
  public:
    RightImpedanceMomentumEq(const Tube&,TubeVertex&,vd& Z);
    ~RightImpedanceMomentumEq();
    vd Error();
    dmat dUi();
    dmat dUim1();
    vd& Z;
    
  }; 

  class RightImpedance:public TubeVertex // Adiabatic impedance boundary condition
  {
  public:
    RightImpedance(const Tube& t,vd Z);
    ~RightImpedance();
    vd Z;			// The impedance
    RightImpedanceMomentumEq mright; // Completely adjusted equation

  };



} // namespace tube

