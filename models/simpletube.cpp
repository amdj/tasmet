#include "bc.h"
#include "models.h"
#include "enginesystem.h"
#include "isentropictube.h"
#include "hopkinslaminarduct.h"
SPOILNAMESPACE
using namespace tasystem;
using namespace tube;

Solver* SimpleTube(us gp,us Nf,d freq,d L,d r,d Tl,d Tr,vd p1,int loglevel,d kappa,int options){
  d S=1;
  cout << "Simpletube called. Options are:\n   ISENTROPIC\n   BLAPPROX\n   DRIVEN\n";
  cout << "Chosen options:\n";
  if(options & ISENTROPIC)
    cout << "ISENTROPIC\n";
  if(options & BLAPPROX)
    cout << "BLAPPROX\n";
  if(options & DRIVEN)
    cout << "DRIVEN\n";
  
  inittrace(loglevel);
  Globalconf air=Globalconf::airSTP(Nf,freq);
  Geom geom1;
  if(options & BLAPPROX)
    geom1=Geom::CylinderBlApprox(gp,L,r);
  else
    geom1=Geom::Cylinder(gp,L,r);
  
  var pL(air,0);
  pL.set(p1);
  tube::LeftPressure bcleft(pL);
  tube::RightIsoTWall bcright(Tr);
  tube::RightAdiabaticWall raw;  
  Tube* t1;
  if(options & ISENTROPIC)
    t1=new IsentropicTube(geom1);
  else
    t1=new HopkinsLaminarDuct(geom1,Tl,Tr);
  
  if(options & DRIVEN)
    t1->addBc(bcleft);
  if(options & ISOTWALL)
    t1->addBc(bcright);
  else
    t1->addBc(raw);

  
  if(options & DRIVEN){
    TaSystem sys(air);
    sys.addSeg(*t1); 
    return new Solver(sys);
  } else {
    EngineSystem sys(air);
    sys.addSeg(*t1); 
    return new Solver(sys);
  }
  delete t1;
}
