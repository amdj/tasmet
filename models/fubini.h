#pragma once
#ifndef _FUBINI_H_
#define _FUBINI_H_

#include "solver.h"
#include "gas.h"

using namespace segment;
using namespace tasystem;
using namespace tube;
using namespace variable;
using namespace gases;
Solver* Fubini(us gp,us Nf,d freq,d L,d S,c p1,int loglevel,d kappa)
{
  initlog(loglevel);
  d dx=L/gp;
  d T0=293.15;
  d p0=101325;
  Gas g("air");
  d rho0=g.rho(T0,p0);
  d c0=g.cm(T0);
  d Z0=rho0*c0;
  d M=pow(abs(p1),2)/p0;
  d phi=1.0;
  d R=sqrt(S/number_pi);
  d PI=S/(2*number_pi*R);
  d rh=S/PI;
  d Mass=0;
  us Ns=2*Nf+1;
  Globalconf gc(Nf,freq,"air",T0,p0,M,S,dx,Mass,kappa);
  Geom geom1(gp,L,S,phi,rh,"inviscid");
  Tube t1(geom1);
  var pL(gc);
  if(Nf>0){
    pL.set(real(p1),1);
    pL.set(imag(p1),2);
  }
  LeftPressure pleft(0,pL);
  
  TAsystem sys(gc);
  sys.addseg(t1);
  sys.addbc(pleft);

  vd Z=Z0*vd(Ns,fillwith::ones);
  RightImpedance iright(0,Z);
  sys.addbc(iright);
  Solver* Sol=new Solver(sys);
  return Sol;  
}



#endif /* _FUBINI_H_ */
