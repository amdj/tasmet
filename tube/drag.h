#pragma once
#ifndef _DRAG_H_
#define _DRAG_H_

#include "var.h"
#include <math_common.h>
#include "globalconf.h"
#include "rottfuncs.h"

//Drag parameter (flow resistance in a tube)

namespace tube{
  SPOILNAMESPACE
  class Tube;
  class TubeVertex;
  class Vertex;
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
    virtual vd operator()(const Vertex&) const;
    virtual dmat dUi(const Vertex&) const;
  protected:
    const Tube& tube;
    const tasystem::Globalconf* gc=NULL;

  };

  class LaminarDragResistance:public DragResistance
  {
  public:
    LaminarDragResistance(const Tube& t);
    vd operator()(const TubeVertex& vertex) const;
    virtual dmat dUi(const TubeVertex&) const;		// Derivative of drag resistance to volume flow
  private:
    vc ComplexResistancecoef(const TubeVertex&) const; // Returns a complex vector of size Ns with drag resistance coefficients for every nonzero frequency (1..Nf)
    laminardrag::ZerofreqDrag zfd;
    rottfuncs::rottfuncs rf;
  };

} // namespace tube
#endif /* _DRAG_H_ */
