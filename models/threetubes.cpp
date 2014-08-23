#include "models.h"

#include "solver.h"
#include "gas.h"
#include "isentropictube.h"
#include "hopkinslaminarduct.h"
#include "twimpedance.h"
#include "pressurebc.h"
using namespace segment;
using namespace tasystem;
using namespace tube;
using namespace variable;
using namespace gases;

Solver* ThreeTubes(us gp,us Nf,d freq,d L,d S1,d S2,vd p1,int loglevel,d kappa)
{
  initlog(loglevel);
  TRACE(100,"L:"<<L);
  d R1=sqrt(S1/number_pi);
  d R2=sqrt(S2/number_pi);
  Globalconf gc=Globalconf::airSTP(Nf,freq,kappa);
  gc.show();
  // Geom geom1=Geom::Cylinder(gp,L,R1);
  // Geom geom2=Geom::Cylinder(gp,L,R2);

  Geom geom1=Geom::CylinderBlApprox(gp,L,R1);
  Geom geom2=Geom::CylinderBlApprox(gp,L,R2);  
  HopkinsLaminarDuct t1(geom1,gc.T0);
  HopkinsLaminarDuct t2(geom2,gc.T0);  
  HopkinsLaminarDuct t3(geom1,gc.T0);
  // IsentropicTube t1(geom1);
  // IsentropicTube t2(geom2);  
  // IsentropicTube t3(geom1);

  var pL(gc);
  for(us i=0;i<gc.Ns;i++)
    pL.set(i,p1(i));

  LeftPressure pleft(pL);
  TaSystem sys(gc);
  t1.addBc(pleft);
  sys.addSeg(t1);
  sys.addSeg(t2);
  sys.addSeg(t3);  

  sys.connectSegs(0,1,tasystem::SegCoupling::tailhead);
  sys.connectSegs(1,2,tasystem::SegCoupling::tailhead);  

  Solver* Sol=new Solver(sys);
  Sol->sys().show(false);
  Sol->sys().init();
  return Sol;  
}



