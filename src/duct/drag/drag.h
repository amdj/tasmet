#pragma once
#ifndef _DRAG_H_
#define _DRAG_H_
#include <tuple>
#include "vtypes.h"

namespace duct{
  SPOILNAMESPACE
  class Cell;

  namespace drag {
    // Return weight factors for values at left and right vertices to
    // value at cell wall
    std::tuple<d,d> wf(const Cell& v);
    
    // Stub for a DragResistance class, with some auxiliary helper functions
    class DragResistance{
    public:
      virtual vd drag(const Cell& v) const;
      virtual dmat dm(const Cell& v) const;
      // virtual dmat drhoi(const Cell& v) const;
      // virtual dmat dpi(const Cell& v) const;

      // Provide a vector of shear wave numbers for each frequency. The
      // length of this vector is Nf+1, and obviously, the first element
      // of this vector is always zero.
      vd shearWaveNumber(const Cell& v) const;
    };
    
    
  } // namespace drag

} // namespace duct
#endif /* _DRAG_H_ */





















