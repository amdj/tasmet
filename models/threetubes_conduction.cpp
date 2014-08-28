#include "models.h"

#include "solver.h"
#include "gas.h"
#include "isentropictube.h"
#include "hopkinslaminarduct.h"
#include "pressurebc.h"
#include "isotwall.h"

using namespace segment;
using namespace tasystem;
using namespace tube;
using namespace variable;
using namespace gases;

Solver* ThreeTubesConduction(us gp,us Nf,d freq,d L,d S1,d S2,vd p1,int loglevel,d kappa,d Tr)
{
  inittrace(loglevel);
  assert(Tr>0);

  TRACE(100,"L:"<<L);
  d R1=sqrt(S1/number_pi);
  d R2=sqrt(S2/number_pi);

  cout << "S1:" <<S1 << "\n";
  cout << "S2:" <<S1 << "\n";  
  cout << "Tr given:" << Tr << "\n";
  cout << "freq given:" << freq << "\n";
  cout << "kappa given:" << kappa << "\n"; 
  cout << "p1 given:" << p1 << "\n"; 
  cout << "gp given:" << gp << "\n"; 
  Globalconf gc=Globalconf::airSTP(Nf,freq,kappa);
  gc.show();
  // Geom geom1=Geom::Cylinder(gp,L,R1);
  // Geom geom2=Geom::Cylinder(gp,L,R2);

  Geom geom1=Geom::CylinderBlApprox(gp,L,R1);
  Geom geom2=Geom::Cylinder(gp,L,R2);  
  // TRACE(100,"ThreetubesConduction hacked to contant cs-area case");
  HopkinsLaminarDuct t1(geom1,gc.T0);
  HopkinsLaminarDuct t2(geom2,gc.T0,Tr);  
  HopkinsLaminarDuct t3(geom1,Tr);

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
  return Sol;  
}



