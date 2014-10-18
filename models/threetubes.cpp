#include "models.h"
#include "geomhelpers.h"
#include "solver.h"
#include "tasystem.h"
#include "gas.h"
#include "isentropictube.h"
#include "hopkinslaminarduct.h"
#include "bc.h"
#include "grid.h"
using namespace segment;
using namespace tasystem;
using namespace tube;
using namespace variable;
using namespace gases;

Solver* ThreeTubes(us gp,us Nf,d freq,d p0,d L,d R1,d R2,vd p1,int loglevel,d kappa,d Tr,int options)
{
  clearConsole();
  inittrace(loglevel);
  TRACE(100,"L:"<<L);
  Globalconf gc=Globalconf::airSTP(Nf,freq,kappa);
  gc.p0=p0;
  gc.show();

  d Sf1=number_pi*pow(R1,2);
  d Sf2=number_pi*pow(R2,2);

  Grid g1(gp,L);
  if(options & BLAYER){
    g1.setLeftBl(5e-6,1.5,40);
    g1.setRightBl(5e-6,1.5,40);
  }
    
  
  Geom geom1,geom2,geom3;
  if(options&BLAPPROX)  {
    geom1=Geom::CylinderBlApprox(g1,R1);
    geom2=Geom::CylinderBlApprox(g1,R2);  
    geom3=Geom::CylinderBlApprox(g1,R1);
  }
  else{
    geom1=Geom::Cylinder(g1,R1);
    geom2=Geom::Cylinder(g1,R2);  
    geom3=Geom::Cylinder(g1,R1);
  }
  segment::smoothEnds(geom1,LAST,geom2,FIRST,5);
  segment::smoothEnds(geom3,FIRST,geom2,LAST,5);
  IsentropicTube t1is(geom1);
  IsentropicTube t2is(geom2);  
  IsentropicTube t3is(geom3);

  HopkinsLaminarDuct t1(geom1,gc.T0);
  HopkinsLaminarDuct t2(geom2,gc.T0,Tr);  
  HopkinsLaminarDuct t3(geom3,Tr);

  var pL(gc);
  for(us i=0;i<gc.Ns();i++)
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
    
  cout << "Tr:" << Tr<<"\n";
  TubeBcVertex* b;
  if(options & ISOTWALL)
    cout << "RightIsoTWall\n";
  if(options & ISOTWALL)
    b=new RightIsoTWall(Tr);
  else
    b=new RightAdiabaticWall();
  t3.addBc(*b);
  t3is.addBc(*b);
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
