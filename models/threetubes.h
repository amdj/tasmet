#pragma once
#ifndef _THREETUBES_H_
#define _THREETUBES_H_

#include "solver.h"
#include "gas.h"
#include "twimpedance.h"
#include "pressurebc.h"
using namespace segment;
using namespace tasystem;
using namespace tube;
using namespace variable;
using namespace gases;
Solver* ThreeTubes(us gp,us Nf,d freq,d L,d S,vd p1,int loglevel,d kappa)
{
  initlog(loglevel);
  d dx=L/gp;
  d T0=293.15;
  d p0=101325;
  Gas g("air");
  d rho0=g.rho(T0,p0);
  d c0=g.cm(T0);
  d z0=rho0*c0;
  d M=pow(max(abs(p1)),2)/p0;
  d phi=1.0;
  d R=sqrt(S/number_pi);
  d Mass=0;
  us Ns=2*Nf+1;
  cout << "Kappa: " << kappa << "\n";
  Globalconf gc(Nf,freq,"air",T0,p0,M,S,dx,Mass,kappa);
  gc.show();
  Geom geom1=Geom::Cylinder(gp,L,R);
  Tube t1(geom1);
  
  // TRACE(30,"p1:"<<p1);
  var pL(gc);
  for(us i=0;i<Ns;i++)
    pL.set(i,p1(i));

  LeftPressure pleft(0,pL);
  TAsystem sys(gc);
  sys.addseg(t1);
  sys.addseg(t1);
  sys.connectSegs(0,1,tasystem::SegCoupling::tailhead);
  sys.addbc(pleft);
  // vd Z=(z0/S)*vd(Ns,fillwith::ones);
  // RightImpedance iright(0,Z);
  // TwImpedance iright(1);
  // sys.addbc(iright);
  Solver* Sol=new Solver(sys);
  // Sol->sys->show();
  return Sol;  
}



#endif /* _THREETUBES_H_ */
