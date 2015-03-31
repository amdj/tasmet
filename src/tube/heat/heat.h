#pragma once
#ifndef _HEAT_H_
#define _HEAT_H_
#include "vtypes.h"

namespace tube{
  SPOILNAMESPACE
  class Cell;
  // Stub for a DragResistance class
  class HeatSource{
  public:
    virtual vd heat(const Cell& v) const;
    virtual dmat dUi(const Cell& v) const;
    virtual dmat dTi(const Cell& v) const;
    // virtual dmat drhoi(const Cell& v) const;
    // virtual dmat dpi(const Cell& v) const;

    // See if comment automatically line on new line das dsag dsag
    // dsag gdsa
    
  };

} // namespace tube
#endif /* _HEAT_H_ */
