#pragma once
#ifndef _DRAG_H_
#define _DRAG_H_

#include "segbase.h"
#include "tubevertex.h"
#include "rottfuncs.h"

//Drag parameter (flow resistance in a tube)

namespace tube{
  SPOILNAMESPACE
  using segment::SegBase;
  using segment::Vertex;

  
  namespace laminardrag{
    // Resistance force for laminar flow for the zero-frequency. 
    d zerodrag_vert(d mu,d rh);
    d zerodrag_circ(d mu,d rh);
    d zerodrag_blapprox(d mu,d rh);
    class ZerofreqDrag{
    public:
      ZerofreqDrag(const SegBase& t);
      d operator()(d mu,d rh,d U) const ; // Full resistance in Newtons
      d operator()(d mu,d rh) const;	// Resistance coefficient
      ~ZerofreqDrag();
    private:
      d (*zerodrag_funptr)(d,d);
      const SegBase& tube;
    };
  } // namespace laminardrag


  class DragResistance{
  public:
    DragResistance(const SegBase& t);
    virtual ~DragResistance();
    virtual vd operator()(const Vertex&) const;
    virtual dmat dUi(const Vertex&) const;
  protected:
    const SegBase& tube;
    const tasystem::Globalconf* gc=NULL;

  };

  class LaminarDragResistance:public DragResistance
  {
  public:
    LaminarDragResistance(const SegBase& t);
    vd operator()(const Vertex& vertex) const;
    virtual dmat dUi(const Vertex&) const;		// Derivative of drag resistance to volume flow
  private:
    vc ComplexResistancecoef(const Vertex&) const; // Returns a complex vector of size Ns with drag resistance coefficients for every nonzero frequency (1..Nf)
    laminardrag::ZerofreqDrag zfd;
    rottfuncs::rottfuncs rf;
  };

} // namespace tube
#endif /* _DRAG_H_ */
