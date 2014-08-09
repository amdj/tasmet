#pragma once
#ifndef _CONETUBE_H_
#define _CONETUBE_H_

#include "solver.h"
#include "gas.h"
#include "hopkinslaminarduct.h"
#include "isentropictube.h"
#include "twimpedance.h"
#include "pressurebc.h"
using namespace segment;
using namespace tasystem;
using namespace tube;
using namespace variable;
using namespace gases;
Solver* ConeTube(us gp,us Nf,d freq,d L,d r1,d r2,vd p1,int loglevel,d kappa)
{
  initlog(loglevel);
  d dx=L/gp;
  d T0=293.15;
  d p0=101325;
  d phi=1.0;

  cout << "Kappa: " << kappa << "\n";
  d S=number_pi*pow(r1,2);
  Globalconf gc(Nf,freq,"air",T0,p0,0,kappa);
  gc.show();
  
  Geom geom1(Geom::Cone(gp,L,r1,r2));
  // HopkinsLaminarDuct t1=HopkinsLaminarDuctTs(geom1,T0);
  IsentropicTube t1(geom1);
  TRACE(100,"Set to isentropic. Still problem with HopkinsLaminarDuct for temperature stuff.");
  // TRACE(30,"p1:"<<p1);
  
  var pL(gc);
  for(us i=0;i<gc.Ns;i++)
    pL.set(i,p1(i));

  LeftPressure pleft(0,pL);
  TwImpedance rightbc(0);
  t1.addBc(pleft);

  TAsystem sys(gc);
  sys.addSeg(t1);


  // sys.addBc(rightbc);
  // vd Z=(z0/S)*vd(Ns,fillwith::ones);
  // RightImpedance iright(0,Z);
  // TwImpedance iright(1);
  // sys.addbc(iright);
  sys.show();
  Solver* Sol=new Solver(sys);
  // Sol->sys->show();
  return Sol;  
}



#endif /* _THREETUBES_H_ */
