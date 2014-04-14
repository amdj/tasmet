// file: bcvertex.h, created March 20th, 2014
// Author: J.A. de Jong

// bcvertex.h: external boundary conditions for tubes. This file
// contains the implementation of typical external boundary conditions
// for tubes as a custom vertex. Examples are adiabatic walls, isothermal walls and an
// adiabatic open pressure boundary conditions.

#include "tube.h"
#include "../var/var.h"
#include "vertex.h"
#include "momentumeq.h"
namespace tube{
  class TubeBcVertex:public TubeVertex
  {
  public:
    TubeBcVertex(const Tube&t,us vertexnr);
    virtual ~TubeBcVertex();
  };

  
  class LeftPressure:public TubeBcVertex
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
    RightImpedanceMomentumEq(const Tube&,TubeBcVertex&,vd& Z);
    ~RightImpedanceMomentumEq();
    vd Error();
    dmat dUi();
    dmat dUim1();
    dmat dpi();
    dmat dpim1();
    dmat drhoim1();
    dmat drhoi();
    vd& Z;
    
  }; 

  class RightImpedance:public TubeBcVertex // Adiabatic impedance boundary condition
  {
  public:
    RightImpedance(const Tube& t,vd Z);
    ~RightImpedance();
    vd Z;			// The impedance
    RightImpedanceMomentumEq mright; // Completely adjusted equation

  };



} // namespace tube

