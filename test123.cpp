/* test.cpp */

#include "tube.h"
#include "conetube.h"
#include "isentropictube.h"
#include "hopkinslaminarduct.h"
#include "var.h"
// #include "enginesystem.h"
#include "tasystem.h"
#include "solver.h"
#include "grid.h"
#include <limits>
#include <boost/archive/text_oarchive.hpp>

using namespace std;
using namespace segment;
using namespace tasystem;
using namespace tube;
TRACETHIS;

int main(int argc,char* argv[]) {
  // clearConsole();
  cout <<  "Running test..." << endl;
  int loglevel=1;
  us gp=4;
  us Nf=0;
  us Ns=2*Nf+1;
  cout << "Max value of int: " << std::numeric_limits<int>::max() << "\n";
  double f=100;
  double omg=2*number_pi*f;
  double T=1/f;
  cout << "SFSG\n";


  if(argc>1)
    loglevel=atoi(argv[1]);
  if(argc>2)
    Nf=atoi(argv[2]);
  if(argc>3)
    gp=atoi(argv[3]);
  inittrace(loglevel);

  cout<< "Loglevel:"<<loglevel<<"\n";

  d L=1;
  VARTRACE(15,L);
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

  // HopkinsLaminarDuct t1(geom1,gc.T0,gc.T0);

  IsentropicTube t1(geom1);
  variable::var pL(gc,0);
  pL.set(0,3.14);
  if(Nf>0)
    pL.set(1,1.0);
  // t1.init(gc);
  // t1.show(1);
  // tube::LeftPressure bcleft(pL);
  // tube::LeftAdiabaticWall bcleft;
  // tube::RightAdiabaticWall bcright;
  // tube::TwImpedance bcright;  
  // tube::RightIsoTWall bcright(T0);
  // TRACE(100,"Add bc to tube...");
  // t1.addBc(bcleft);
  // t1.addBc(bcright);


  // EngineSystem sys(air);
  // sys.setTimingConstraint(0,0,3,2);
  // sys.setAmplitudeDof(0,0,3,1);
  TaSystem sys(air);
  sys.addSeg(t1); 
  TRACE(25,"SFSG");
  static_cast<Tube*>(sys.getSeg(0))->geom().show();
  sys.init();

  // Solver sol(sys);

  sys.show(2);
  // sol.sys().showJac();

  // // vd domg(15,fillwith::zeros);
  // // sol.sys().getSeg(0)->domg(domg);
  // // cout << "domg: "<< domg << "\n";
  // cout << "Result:\n"<< sol.sys().getRes() << "\n";
  // cout << "Error:\n"<< sol.sys().error() << "\n";

  // sol.doIter();
  // static_cast<Tube*>(sol.sys().getSeg(0))->getResAt(0,0);
  return 0;
}










