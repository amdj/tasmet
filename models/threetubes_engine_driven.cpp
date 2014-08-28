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

inline us max(us x,us y){ return x>y?x:y;}

Solver* ThreeTubesEngineDriven(us gp,us Nf,d freq,d Tr,vd p1,int loglevel,d kappa)
{
  TRACE(30,"Entering ThreeTubesEngine");
  inittrace(loglevel);
  assert(Tr>0);

  d p0=376e3;
  d T0=293.15;


  d L1=90e-2;
  d L2=35e-3;
  d L3=1-L1-L2;

  d Ltot=L1+L2+L3;
  
  us gp1=max(round(gp*L1/Ltot),4);  
  us gp2=max(round(gp*L2/Ltot),4)*5;
  us gp3=max(round(gp*L3/Ltot),4)*2;
  
  cout << "gp1: "<< gp1<< "\n";
  cout << "gp2: "<< gp2<< "\n";
  cout << "gp3: "<< gp3<< "\n";
  d Rtube=38.2e-3/2;
  d S=number_pi*pow(Rtube,2);
  
  Globalconf gc(Nf,freq,"helium",T0,p0,kappa);
  gc.show();

  d phi_s=0.73;
  d y0=0.77e-3/2;
  
  Geom geom1=Geom::CylinderBlApprox(gp1,L1,Rtube);
  Geom geom2=Geom::VertPlates(gp2,L2,S,phi_s,y0);
  Geom geom3=Geom::CylinderBlApprox(gp3,L3,Rtube);

  HopkinsLaminarDuct t1(geom1,gc.T0);
  HopkinsLaminarDuct t2(geom2,gc.T0,Tr);  
  HopkinsLaminarDuct t3(geom3,Tr);

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
  return Sol;  
}



