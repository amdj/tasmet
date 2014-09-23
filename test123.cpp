/* test.cpp */

#include "tube.h"
#include "isentropictube.h"
#include "hopkinslaminarduct.h"
#include "enginesystem.h"
#include "solver.h"
#include "bc.h"
#include <limits>
using namespace std;
using namespace segment; 
using namespace tasystem;
using namespace tube;


int main(int argc,char* argv[]) {
  cout <<  "Running test..." << endl;
  int loglevel=25;
  us gp=4;
  us Nf=0;
  us Ns=2*Nf+1;

  cout << "Max value of int: " << std::numeric_limits<int>::max() << "\n";
  double f=100;
  double omg=2*number_pi*f;
  double T=1/f;
  if(argc>1)
    loglevel=atoi(argv[1]);
  if(argc>2)
    Nf=atoi(argv[2]);
  if(argc>3)
    gp=atoi(argv[3]);
  cout<< "Loglevel:"<<loglevel<<"\n";
  inittrace(loglevel);
  d L=0.01;

  d S=1;
  
  d rtube=sqrt(S/number_pi);

  d phi=1.0;
  d PI=2*number_pi*rtube;
  d rh=S/PI;
  d p0=101325.0;
  d T0=293.15;
  d griddx=L/gp;
  d S0=S;


  d kappa=0.1;
  Globalconf air=Globalconf::airSTP(Nf,f);
  Globalconf gc=air;
  Geom geom1=Geom::CylinderBlApprox(gp,L,rtube);
  HopkinsLaminarDuct t1(geom1,gc.T0,gc.T0);
  // IsentropicTube t1(geom1);
  var pL(gc,0);
  // pL.set(0,3.14);
  if(Nf>0)
    pL.set(1,1.0);
  // tube::LeftPressure bcleft(pL);
  tube::LeftAdiabaticWall bcleft;
  // tube::RightAdiabaticWall bcright;
  // tube::TwImpedance raw;  
  tube::RightIsoTWall bcright(T0+10);
  // TRACE(100,"Add bc to tube...");
  t1.addBc(bcleft);
  t1.addBc(bcright);


  // EngineSystem sys(air);
  // sys.setTimingConstraint(0,0,3,2);
  // sys.setAmplitudeDof(0,0,3,1);
  EngineSystem sys(air);
  // TaSystem sys(air);
  sys.addSeg(t1); 
  sys.init();

  Solver sol(sys);
  sol.sys().show(true);
  sol.sys().showJac();

  // vd domg(15,fillwith::zeros);
  // sol.sys().getSeg(0)->domg(domg);
  // cout << "domg: "<< domg << "\n";
  cout << "Result:\n"<< sol.sys().getRes() << "\n";
  cout << "Error:\n"<< sol.sys().error() << "\n";

  sol.doIter();
  return 0;
}

