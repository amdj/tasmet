#pragma once
#ifndef _PICKADOF_H_
#define _PICKADOF_H_
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
      for(us i=0;i<vertexnr;i++)
        dofnr+=static_cast<segment::Seg*>(sys.getSeg(segnr))->vvertex.at(i)->getNDofs();
      if(varnr>0)
        dofnr+=(varnr-1)*sys.gc.Ns;
      dofnr+=freqnr;
      return dofnr;
    }
  };

}		// namespace tasystem
#endif /* _PICKADOF_H_ */
