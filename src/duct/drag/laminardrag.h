#include "drag.h"
#include "rottfuncs.h"

namespace duct{
  SPOILNAMESPACE
  namespace drag {
    
    // The Drag coefficient for "frequency" zero.
    class ZeroFreqDragCoef;    

    // Laminar drag resistance
    class LaminarDragResistance:public DragResistance
    {
      ZeroFreqDragCoef* zfd;
      rottfuncs::RottFuncs rf;
    public:
      LaminarDragResistance(const Duct& t);
      LaminarDragResistance(const LaminarDragResistance&)=delete;
      LaminarDragResistance& operator=(const LaminarDragResistance&)=delete;
      ~LaminarDragResistance();

      // Overloaded virtuals
      vd drag(const Cell& cell) const;
      tasystem::JacRow dDrag(const Cell&) const;	       

    private:
      dmat dm(const Cell&) const;		// Derivative of drag resistance
      // to volume flow
    
      // Returns a complex vector of size Ns with drag resistance
      // coefficients for every nonzero frequency (1..Nf)
      vc ComplexResistancecoef(const Cell&) const;

    };
  } // namespace drag
} // namespace duct
