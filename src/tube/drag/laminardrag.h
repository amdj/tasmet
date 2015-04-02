#include "drag.h"
#include "rottfuncs.h"

namespace tube{
  SPOILNAMESPACE
  namespace  {
    class ZeroFreqDrag;    
  } // namespace 

  
  class LaminarDragResistance:public DragResistance
  {
    ZeroFreqDrag* zfd;
    rottfuncs::RottFuncs rf;
  public:
    LaminarDragResistance(const Tube& t);
    LaminarDragResistance(const LaminarDragResistance&)=delete;
    LaminarDragResistance& operator=(const LaminarDragResistance&)=delete;
    ~LaminarDragResistance();

    // Overloaded virtuals
    vd drag(const Cell& cell) const;
    dmat dm(const Cell&) const;		// Derivative of drag resistance
                                    // to volume flow
    
  private:
    // Returns a complex vector of size Ns with drag resistance
    // coefficients for every nonzero frequency (1..Nf)
    vc ComplexResistancecoef(const Cell&) const;

  };
} // namespace tube
