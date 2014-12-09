#pragma once
#ifndef _DRAG_H_
#define _DRAG_H_
#include "vtypes.h"

namespace tube{
  SPOILNAMESPACE
  class TubeVertex;
  // Stub for a DragResistance class
  class DragResistance{
  public:
    virtual vd drag(const TubeVertex& v) const;
    virtual dmat dUi(const TubeVertex& v) const;
    // virtual dmat drhoi(const TubeVertex& v) const;
    // virtual dmat dpi(const TubeVertex& v) const;

  };

} // namespace tube
#endif /* _DRAG_H_ */





















