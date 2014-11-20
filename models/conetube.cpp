#include "solver.h"
#include "tasystem.h"
#include "gas.h"
#include "conetube.h"
#include "hopkinslaminarduct.h"
#include "isentropictube.h"
#include "twimpedance.h"
#include "pressurebc.h"
using namespace segment;
using namespace geom;
using namespace tasystem;
using namespace tube;
using namespace variable;
using namespace gases;
Solver* ConeTubeModel(us gp,us Nf,d freq,d L,d r1,d r2,vd p1,int loglevel,d kappa,us isentropic,us blapprox)
{
  inittrace(loglevel);
  d phi=1.0;
  cout << "Kappa: " << kappa << "\n";
  d S=number_pi*pow(r1,2);
  Globalconf gc=Globalconf::airSTP(Nf,freq);
  Geom* geom1;
  Grid grid(gp,L);
  if(blapprox)
    geom1=new ConeTube(grid,r1,r2,true);
  else
    geom1=new ConeTube(grid,r1,r2,false);
    
  TaSystem sys(gc);
  var pL(gc);
  for(us i=0;i<gc.Ns();i++)
    pL.set(i,p1(i));
  LeftPressure pleft(pL);
  // TwImpedance rightbc;


  if(isentropic){
    IsentropicTube t1(geom1);
    t1.addBc(pleft);
    // t1.addBc(rightbc);
    sys.addSeg(t1);
  }
  else{
    HopkinsLaminarDuct t1(geom1,gc.T0);
    t1.addBc(pleft);
    // t1.addBc(rightbc);    
    sys.addSeg(t1);
  }

  Solver* Sol=new Solver(sys);
  return Sol;  
}

