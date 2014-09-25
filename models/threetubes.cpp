#include "models.h"

#include "solver.h"
#include "tasystem.h"
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
  gc.setp0(50*gc.p0);
  gc.show();

  d Sf1=number_pi*pow(R1,2);
  d Sf2=number_pi*pow(R2,2);
  Geom geom1,geom2,geom3;
  if(options&BLAPPROX)  {
    geom1=Geom::CylinderBlApprox(gp,L,R1);
    geom2=Geom::CylinderBlApprox(gp,L,R2);  
    geom3=Geom::CylinderBlApprox(gp,L,R1);
  }
  else{
    geom1=Geom::Cylinder(gp,L,R1);
    geom2=Geom::Cylinder(gp,L,R2);  
    geom3=Geom::Cylinder(gp,L,R1);
  }
  segment::smoothEnds(geom1,LAST,geom2,FIRST);
  segment::smoothEnds(geom2,LAST,geom3,FIRST);
  IsentropicTube t1is(geom1);
  IsentropicTube t2is(geom2);  
  IsentropicTube t3is(geom3);

  HopkinsLaminarDuct t1(geom1,gc.T0);
  HopkinsLaminarDuct t2(geom2,gc.T0,Tr);  
  HopkinsLaminarDuct t3(geom3,Tr);

  var pL(gc);
  for(us i=0;i<gc.Ns;i++)
    pL.set(i,p1(i));
  if(options&DRIVEN){
    cout << "Driven system\n";
    LeftPressure pleft(pL);
    t1.addBc(pleft);
    t1is.addBc(pleft);
  }
  else{
    cout << "Non-driven system\n";
    LeftIsoTWall pleft(gc.T0);
    t1.addBc(pleft);
    t1is.addBc(pleft);
  }
    
  // RightAdiabaticWall right;
  cout << "Tr:" << Tr<<"\n";
  TubeBcVertex* b;
  RightAdiabaticWall raw;
  if(options & BLAPPROX)
    b=new RightAdiabaticWall();
  else
    b=new RightIsoTWall(Tr);
  t3.addBc(raw);
  t3is.addBc(raw);
  delete b;
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
