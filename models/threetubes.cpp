#include "models.h"

#include "solver.h"
#include "gas.h"
#include "isentropictube.h"
#include "hopkinslaminarduct.h"
#include "bc.h"
using namespace segment;
using namespace tasystem;
using namespace tube;
using namespace variable;
using namespace gases;

Solver* ThreeTubes(us gp,us Nf,d freq,d L,d R1,d R2,vd p1,int loglevel,d kappa,d Tr,int options)
{
  inittrace(loglevel);
  TRACE(100,"L:"<<L);
  Globalconf gc=Globalconf::airSTP(Nf,freq,kappa);
  gc.show();

  Geom geom1=Geom::CylinderBlApprox(gp,L,R1);
  Geom geom2=Geom::CylinderBlApprox(gp,L,R2);  
  IsentropicTube t1is(geom1);
  IsentropicTube t2is(geom2);  
  IsentropicTube t3is(geom1);

  HopkinsLaminarDuct t1(geom1,gc.T0);
  HopkinsLaminarDuct t2(geom2,gc.T0);  
  HopkinsLaminarDuct t3(geom1,gc.T0);

  var pL(gc);
  for(us i=0;i<gc.Ns;i++)
    pL.set(i,p1(i));

  LeftPressure pleft(pL);
  t1.addBc(pleft);
  t1is.addBc(pleft);
  
  RightAdiabaticWall right;
  // RightIsoTWall right(Tr);
  t3.addBc(right);
  t3is.addBc(right);
  TaSystem sys(gc);
  if(options & ISENTROPIC){
    sys.addSeg(t1is);
    sys.addSeg(t2is);
    sys.addSeg(t3is);  
  }
  else{
    sys.addSeg(t1);
    sys.addSeg(t2);
    sys.addSeg(t3);  
  }
  sys.connectSegs(0,1,tasystem::SegCoupling::tailhead);
  sys.connectSegs(1,2,tasystem::SegCoupling::tailhead);  

  Solver* Sol=new Solver(sys);
  return Sol;  
}
