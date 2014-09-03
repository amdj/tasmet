#include "bc.h"
#include "solver.h"
#include "enginesystem.h"
#include "tube.h"
#include "hopkinslaminarduct.h"
SPOILNAMESPACE
using namespace tasystem;
using namespace tube;

Solver* SimpleTube(us gp,us Nf,d freq,d L,d r,d Tl,d Tr,vd p1,int loglevel,d kappa,us blapprox,us driven,us rwall){
  d S=1;
  inittrace(loglevel);
  Globalconf air=Globalconf::airSTP(Nf,freq);
  Geom geom1;
  if(blapprox)
    geom1=Geom::CylinderBlApprox(gp,L,r);
  else
    geom1=Geom::Cylinder(gp,L,r);
  
  var pL(air,0);
  pL.set(p1);
  tube::LeftPressure bcleft(pL);
  tube::RightIsoTWall bcright(Tr);
  
  HopkinsLaminarDuct t1(geom1,Tl,Tr);
  if(driven)
    t1.addBc(bcleft);
  if(rwall)
    t1.addBc(bcright);

  if(driven){
    TaSystem sys(air);
    sys.addSeg(t1); 
    return new Solver(sys);
  } else {
    EngineSystem sys(air);
    sys.addSeg(t1); 
    return new Solver(sys);
  }
  
}
