#pragma once
#ifndef _HEAT_H_
#define _HEAT_H_
#include <vtypes.h>

namespace tube{
  SPOILNAMESPACE
  class TubeVertex;
  // Stub for a DragResistance class
  class HeatSource{
  public:
    virtual vd heat(const TubeVertex& v) const;
    virtual dmat dUi(const TubeVertex& v) const;
    virtual dmat dTi(const TubeVertex& v) const;
    // virtual dmat drhoi(const TubeVertex& v) const;
    // virtual dmat dpi(const TubeVertex& v) const;

    // See if comment automatically line on new line das dsag dsag
    // dsag gdsa
    
  };

} // namespace tube
#endif /* _HEAT_H_ */
