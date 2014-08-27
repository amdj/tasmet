#pragma once
#ifndef _TIMINGCONSTRAINT_H_
#define _TIMINGCONSTRAINT_H_
#include "vtypes.h"
#include "tasystem.h"
#include "seg.h"
namespace tasystem{
  SPOILNAMESPACE
  class PickADof{
    us segnr,vertexnr,varnr,freqnr;
  public:
    PickADof(us segnr=0,us vertexnr=0,us varnr=0,us freqnr=0)
    {set(segnr,vertexnr,varnr,freqnr);}
    void set(us segnr,us vertexnr,us varnr,us freqnr){
      this->segnr=segnr;
      this->vertexnr=vertexnr;
      this->varnr=varnr;
      this->freqnr=freqnr;
    }
    d value(const TaSystem& sys) const {
      TRACE(15,"PickADof::error()");
      us Ns=sys.gc.Ns;
      return static_cast<segment::Seg*>(sys.getSeg(segnr))->vvertex.at(vertexnr)->getRes()(varnr*Ns+freqnr);
    }
    us dofnr(const TaSystem& sys) const {
      TRACE(15,"PickADof::dofnr()");
      us dofnr=0;
      for(us i=0;i<segnr;i++)
	dofnr+=sys.getSeg(i)->getNDofs();
      dofnr+=vertexnr*Neq*sys.gc.Ns+varnr*sys.gc.Ns+freqnr;
      return dofnr;
    }
  };

}		// namespace tasystem
#endif /* _TIMINGCONSTRAINT_H_ */
