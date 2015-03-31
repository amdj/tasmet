#pragma once
#ifndef _DRAG_H_
#define _DRAG_H_
#include "vtypes.h"

namespace tube{
  SPOILNAMESPACE
  class Cell;
  // Stub for a DragResistance class
  class DragResistance{
  public:
    virtual vd drag(const Cell& v) const;
    virtual dmat dUi(const Cell& v) const;
    // virtual dmat drhoi(const Cell& v) const;
    // virtual dmat dpi(const Cell& v) const;

  };

} // namespace tube
#endif /* _DRAG_H_ */





















