/* test.cpp */

#include "tube.h"
#include "conetube.h"
#include "isentropictube.h"
#include "hopkinslaminarduct.h"
#include "var.h"

#include "adiabaticwall.h"
#include "isotwall.h"
#include "pressurebc.h"
// #include "enginesystem.h"

#include "tasystem.h"
#include "solver.h"
#include "grid.h"


SPOILNAMESPACE
using namespace segment;
using namespace tasystem;
using namespace tube;


int main(int argc,char* argv[]) {
  // clearConsole();
  cout <<  "Running test..." << endl;
  int loglevel=20;
  us gp=4;
  us Nf=0;
  us Ns=2*Nf+1;
  // cout << "Max value of int: " << std::numeric_limits<int>::max() << "\n";
  double f=100;
  double omg=2*number_pi*f;
  double T=1/f;

  if(argc>1)
    loglevel=atoi(argv[1]);
  if(argc>2)
    Nf=atoi(argv[2]);
  if(argc>3)
    gp=atoi(argv[3]);
  inittrace(loglevel);

  cout<< "Loglevel:"<<loglevel<<"\n";

  d L=1;
  d S=1;
  d rtube=sqrt(S/number_pi);
  d phi=1.0;
  d PI=2*number_pi*rtube;
  d rh=S/PI;
  d p0=101325.0;
  d T0=293.15;
  d S0=S;


  Grid grid(gp,L);
  // grid.setLeftBl(1e-5,5,40);
  // grid.setRightBl(1e-5,5,40);  
  CylindricalTube geom1(grid,rtube);

  // d kappa=0.1;
  Globalconf air=Globalconf::airSTP(Nf,f);
  // {
  //   std::ofstream ofs("filename.txt");
  //   boost::archive::text_oarchive oa(ofs);
  //   oa << air;
  // }

  Globalconf gc=air;

  HopkinsLaminarDuct t1(geom1,gc.T0(),gc.T0());

  // IsentropicTube t1(geom1);
  tasystem::var pL(gc,0);
  pL.set(0,3.14);
  if(Nf>0)
    pL.set(1,1.0);
  // PressureBc first(pL,0,pos::left);
  // PressureBc p(pL,0,pos::left);
  tasystem::var Tbc(gc,393.15);
  // AdiabaticWall first(0,pos::left);
  PressureBc first(pL,0,left); 
  IsoTWall second(0,right,Tbc);
 
  // AdiabaticWall bright(0,pos::right);
  // EngineSystem sys(air);
  // sys.setTimingConstraint(0,0,3,2);
  // sys.setAmplitudeDof(0,0,3,1);
  TaSystem sys(air);
  sys+=t1;
  sys+=first;
  sys+=second;
  // sys.show(20);

  Solver sol(sys);
  sys.show(10);
  // sys.show(2);
  TRACE(20,"Showing system...");
  // sol.sys().show(10);
  TRACE(20,"Computing error...");

  cout << "Error:\n"<< sys.error() << "\n";
  cout << "Res:\n"<< sys.getRes() << "\n";
  sys.showJac();
  TRACE(20,"Hoi");
  evd error=sys.error();
  esdmat jac=sys.jac();

  // TRACE(20,"Hoi");
  // VARTRACE(15,error);
  // VARTRACE(15,jac);
  // evd solu=tasystem::solvesys_eigen(jac,error);
  sol.doIter();
  // TRACE(20,"Hoi");
  // sleep(1);
  TRACE(15,"err"<<sol.sys().error());
  return 0;
}










