#include "models.h"
#include "solver.h"
#include "gas.h"
#include "bc.h"
#include "isentropictube.h"
#include "hopkinslaminarduct.h"
#include "enginesystem.h"
using namespace segment;
using namespace tasystem;
using namespace tube;
using namespace variable;
using namespace gases;

inline us max(us x,us y){ return x>y?x:y;}


Solver* Atchley_Engine(us gp,us Nf,d freq,d Tr,int loglevel,d kappa,vd p1,d p0,int options)
{
  // #Some global params
  inittrace(loglevel);
  // d p0=376e3;
  d T0=293.15;  
  Globalconf gc(Nf,freq,"helium",T0,p0,kappa);

  d Rtube=38.2e-3/2;
  d S0=number_pi*pow(Rtube,2);
  d rhtube=Rtube/2;
  d Ltot=1.0;
  // ############################## Resonator
  d Lres=87.97e-2;
  us gpres=max(round(gp*Lres/Ltot),4);  
  Geom resgeom=Geom::CylinderBlApprox(gpres,Lres,Rtube);

  // ############################## Cold HX's
  d Lchx=10.2e-3;
  d y0chx=1.02e-3/2;
  d phichx=0.70;
  d Lchxgap=1.5e-3;
  us gphx=25;  
  Geom chxgeom=Geom::VertPlates(gphx,Lchx,S0,phichx,y0chx);

  // ############################## Stack
  d y0stk=0.77e-3/2;
  d Lstk=3.5e-2;
  us gpstk=60;
  d phistk=0.73;
  Geom stkgeom=Geom::VertPlates(gpstk,Lstk,S0,phistk,y0stk);

  // ############################## Hot HX
  d y0hhx=y0chx;
  d Lhhx=7.62e-3;
  d phihhx=0.70;
  Geom hhxgeom=Geom::VertPlates(gphx,Lhhx,S0,phihhx,y0hhx);

  // ############################## Hot end
  d Lhotend=(5e-2)-Lhhx;
  us gphotend=90;
  Geom hotendgeom=Geom::CylinderBlApprox(gphotend,Lhotend,Rtube);

  
  // And blend togeter
  smoothEnds(resgeom,LAST,chxgeom,FIRST);
  smoothEnds(chxgeom,LAST,stkgeom,FIRST);
  smoothEnds(stkgeom,LAST,hhxgeom,FIRST);
  smoothEnds(hhxgeom,LAST,hotendgeom,FIRST);
  
  HopkinsLaminarDuct resonator(resgeom,gc.T0);
  HopkinsLaminarDuct chx(chxgeom,gc.T0);
  HopkinsLaminarDuct stk(stkgeom,gc.T0,Tr);
  HopkinsLaminarDuct hhx(hhxgeom,Tr);
  HopkinsLaminarDuct hotend(hotendgeom,Tr);  

  RightAdiabaticWall rwall;
  hotend.addBc(rwall);
  
  TaSystem sys(gc);
  sys.addSeg(resonator);
  sys.addSeg(chx);
  sys.addSeg(stk);
  sys.addSeg(hhx);  
  sys.addSeg(hotend);  
  sys.connectSegs(0,1,tasystem::SegCoupling::tailhead);
  sys.connectSegs(1,2,tasystem::SegCoupling::tailhead);
  sys.connectSegs(2,3,tasystem::SegCoupling::tailhead);
  sys.connectSegs(3,4,tasystem::SegCoupling::tailhead);  
  
  if(options & DRIVEN){
    cout << "Driven system\n";
    var pL(gc);
    for(us i=0;i<gc.Ns;i++)
      pL.set(i,p1(i));
    LeftPressure l1(pL);
    static_cast<Tube*>(sys.getSeg(0))->addBc(l1);
    
    Solver* Sol=new Solver(sys);
    return Sol;

  }
  else{
    cout << "NOT driven system\n";

    EngineSystem ensys(sys);
    // CHANGED to constraint on second harmonic, see what happens
    // ensys.setTimingConstraint(0,0,3,4);
    ensys.setTimingConstraint(0,0,3,2);
    ensys.setAmplitudeDof(0,0,3,1);
    Solver* Sol=new Solver(ensys);
    return Sol;
  }
}


