#include "drag.h"
#include "rottfuncs.h"
#include "globalconf.h"

namespace tube{
  SPOILNAMESPACE

  
  namespace laminardrag{
    // Resistance force for laminar flow for the zero-frequency. 
    class ZeroFreqDrag{
    public:
      ZeroFreqDrag(const Tube& t);
      d operator()(d mu,d rh,d U) const ; // Full resistance in Newtons
      d operator()(d mu,d rh) const;	// Resistance coefficient
    private:
      d (*zerodrag_funptr)(d,d);
    };
  } // namespace laminardrag



  class LaminarDragResistance:public DragResistance
  {
  public:
    LaminarDragResistance(const Tube& t);
    d fnu(us i) const;
    virtual vd drag(const TubeVertex& vertex) const;
    virtual dmat dUi(const TubeVertex&) const;		// Derivative of drag resistance to volume flow
    vc ComplexResistancecoef(const TubeVertex&) const; // Returns a complex vector of size Ns with drag resistance coefficients for every nonzero frequency (1..Nf)
  private:
    laminardrag::ZeroFreqDrag zfd;
    rottfuncs::RottFuncs rf;
  };
} // namespace tube
