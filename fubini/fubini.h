#pragma once
#ifndef _FUBINI_H_
#define _FUBINI_H_

#include "solver.h"
#include "gas.h"

Solver* Fubini(us gp,us Nf,d freq,d L,d S,c p1,int loglevel,d kappa)
{
  initlog(loglevel);
  d dx=L/gp;
  d T0=293.15;
  d p0=101325;
  d M=pow(abs(p1),2)/p0;

  Globalconf gc(Nf,freq,"air",T0,p0,M,S,dx,Mass,kappa);

  Solver* Sol;
  retur Sol;  
}



#endif /* _FUBINI_H_ */
