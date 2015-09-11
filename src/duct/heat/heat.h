#pragma once
#ifndef _HEAT_H_
#define _HEAT_H_
#include "vtypes.h"

namespace tasystem {
    class JacRow;
} // tasystem

namespace duct{
  SPOILNAMESPACE
  class Cell;
  // Stub for a DragResistance class
  class HeatSource{
  public:
    // Returns heat flow from solid to fluid in W/m of duct length. Should be added
    // possitively to heat balance of fluid and negatively to heat
    // balance of solid
    virtual vd Qsf(const Cell& v) const;
    virtual tasystem::JacRow dQsf(const Cell& v) const;
    
  };

} // namespace duct
#endif /* _HEAT_H_ */
