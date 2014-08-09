#pragma once
#ifndef _FUBINI_H_
#define _FUBINI_H_

#include "solver.h"
#include "gas.h"
#include "twimpedance.h"
#include "isotwall.h"
#include "laminarduct.h"
#include "isentropictube.h"
#include "pressurebc.h"
using namespace segment;
using namespace tasystem;
using namespace tube;
using namespace variable;
using namespace gases;
Solver* Fubini(us gp,us Nf,d freq,d L,d S,vd p1,int loglevel,d kappa)
{
  initlog(loglevel);
  d dx=L/gp;
  d T0=293.15;
  d p0=101325;
  Gas g("air");
  d rho0=g.rho(T0,p0);
  d c0=g.cm(T0);
  d z0=rho0*c0;
  d R=sqrt(S/number_pi);

  d PI=S/(2*number_pi*R);
  d Mass=0;
  us Ns=2*Nf+1;
  cout << "Kappa: " << kappa << "\n";
  Globalconf gc(Nf,freq,"air",T0,p0,Mass,kappa);
  gc.show();
  Geom geom1=Geom::Cylinder(gp,L,R);
  IsentropicTube t1(geom1);
  // TRACE(30,"p1:"<<p1);
  var pL(gc);
  for(us i=0;i<Ns;i++)
    pL.set(i,p1(i));

  vd Zv=(z0/S)*vd(2*Nf+1,fillwith::ones);
  // cout << Zv;
  LeftPressure bleft(0,pL);
  TwImpedance bright(0);
  // RightImpedance bright(0,Zv);
  // RightIsoTWall bright(0,T0);
  TAsystem sys(gc);
  t1.addBc(bleft);
  t1.addBc(bright);

  sys.addSeg(t1);
  Solver* Sol=new Solver(sys);
  Sol->sys.show(false);
  return Sol;  
}



#endif /* _FUBINI_H_ */

















