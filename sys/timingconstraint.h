#pragma once
#ifndef _TIMINGCONSTRAINT_H_
#define _TIMINGCONSTRAINT_H_
#include "vtypes.h"
#include "tasystem.h"
#include "seg.h"
namespace tasystem{
  SPOILNAMESPACE
  class TimingConstraint{
    us segnr,vertexnr,varnr,freqnr;
  public:
    TimingConstraint(us segnr=0,us vertexnr=0,us varnr=3,us freqnr=2):
      segnr(segnr),vertexnr(vertexnr),varnr(varnr),freqnr(freqnr)
    {}
    d error(const TaSystem& sys) const {
      TRACE(15,"TimingConstraint::error()");
      us Ns=sys.gc.Ns;
      return static_cast<segment::Seg*>(sys.getSeg(segnr))->vvertex.at(vertexnr)->getRes()(varnr*Ns+freqnr);
    }
  };

}		// namespace tasystem
#endif /* _TIMINGCONSTRAINT_H_ */
