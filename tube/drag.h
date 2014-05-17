#pragma once
#include "../var/var.h"
#include <math_common.h>
#include "../globalconf.h"

//Drag parameter (flow resistance in a tube)

namespace tube{
  SPOILNAMESPACE
  class Tube;
  namespace laminardrag{
    // Resistance force for laminar flow for the zero-frequency. 
    d zerodrag_vert(d mu,d rh);
    d zerodrag_circ(d mu,d rh);
    d zerodrag_blapprox(d mu,d rh);
    class ZerofreqDrag{
    public:
      ZerofreqDrag(const Tube& t);
      d operator()(d mu,d rh,d U) const ; // Full resistance in Newtons
      d operator()(d mu,d rh) const;	// Resistance coefficient
      ~ZerofreqDrag();
    private:
      d (*zerodrag_funptr)(d,d);
      const Tube& tube;
    };
  } // namespace laminardrag


  class DragResistance{
  public:
    DragResistance(const Tube& t);
    virtual ~DragResistance();
    virtual vd operator()(us i) const;
    virtual dmat dUi(us i) const;
  protected:
    const Tube& tube;
    const tasystem::Globalconf& gc;

  };

  class LaminarDragResistance:public DragResistance
  {
  public:
    LaminarDragResistance(const Tube& t);
    vd operator()(us i) const;
    dmat dUi(us i) const;		// Derivative of drag resistance to volume flow
  private:
    vc ComplexResistancecoef(us i) const; // Returns a complex vector of size Ns with drag resistance coefficients for every nonzero frequency (1..Nf)
    laminardrag::ZerofreqDrag zfd;
    math_common::rottfuncs rf;
  };

} // namespace tube
