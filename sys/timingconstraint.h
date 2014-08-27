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
    us dofnr(const TaSystem& sys) const {
      us dofnr;
      if(segnr==0)
	dofnr=vertexnr*Neq*sys.gc.Ns+varnr*sys.gc.Ns+freqnr;
      else {
	WARN("Extended varnr for not the first segment not yet implmented! Exiting...");
	exit(1);
      }
      return dofnr;
    }
  };

}		// namespace tasystem
#endif /* _TIMINGCONSTRAINT_H_ */
