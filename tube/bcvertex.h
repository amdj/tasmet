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
    LeftPressure(const Tube&t,variable::var& pres,d T0);
    LeftPressure(const Tube&t,variable::var& pres,variable::var& temp);

    ~LeftPressure();
    void Init();
    // virtual vd cbcsource();
    virtual vd msource() const;	// Prescribed left pressure
    virtual vd esource() const;	// Same
    variable::var pL;			// Pressure boundary condition
    variable::var TL;			// Temperature boundary conditions
  };

  class MomentumLeftPressure:public Momentum{
  public:
    MomentumLeftPressure(const Tube& tube,TubeVertex& gp);
    ~MomentumLeftPressure();
    
  };				// MomentumLeftPressure
  




} // namespace tube

