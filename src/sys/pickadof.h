#pragma once
#ifndef _PICKADOF_H_
#define _PICKADOF_H_
#include "vtypes.h"
#include "tasystem.h"
#include "seg.h"
#include "cell.h"

using segment::Seg;

namespace tasystem{
  SPOILNAMESPACE
  class PickADof{
    us segnr,cellnr,varnr,freqnr;
  public:
    PickADof(us segnr=0,us cellnr=0,us varnr=0,us freqnr=0)
    {set(segnr,cellnr,varnr,freqnr);}
    void set(us segnr,us cellnr,us varnr,us freqnr){
      this->segnr=segnr;
      this->cellnr=cellnr;
      this->varnr=varnr;
      this->freqnr=freqnr;
    }
    d value(const TaSystem& sys) const {
      TRACE(15,"PickADof::error()");
      WARN("Not working anymore");
// us Ns=sys.gc.Ns();
      // us i=0;
      // while()
      // Seg* seg=sys.getSeg()[i];
      // return static_cast<segment::Seg*>(sys.getSeg(segnr))->cells.at(cellnr)->getRes()(varnr*Ns+freqnr);
      return 0;
    }
    us dofnr(const TaSystem& sys) const {
      TRACE(15,"PickADof::dofnr()");
      us dofnr=0;
      WARN("Not working anymore");
      // for(us i=0;i<segnr;i++)
      //   dofnr+=sys.getSeg(i)->getNDofs();
      // for(us i=0;i<cellnr;i++)
      //   dofnr+=static_cast<segment::Seg*>(sys.getSeg(segnr))->cells.at(i)->getNDofs();
      // dofnr+=varnr*sys.gc.Ns();
      // dofnr+=freqnr;
      // return dofnr;
      return 0;
    }
  };

}		// namespace tasystem
#endif /* _PICKADOF_H_ */
