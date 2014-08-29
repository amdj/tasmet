/* test.cpp */

#include "tube.h"
#include "isentropictube.h"
#include "hopkinslaminarduct.h"
#include "enginesystem.h"
#include "solver.h"
#include "bc.h"
using namespace std;
using namespace segment; 
using namespace tasystem;
using namespace tube;


int main(int argc,char* argv[]) {
  cout <<  "Running test..." << endl;
  int loglevel=25;
  us gp=4;
  us Nf=1;
  us Ns=2*Nf+1;
  double f=100;
  double omg=2*number_pi*f;
  double T=1/f;
  if(argc>1)
    loglevel=atoi(argv[1]);
  if(argc>2)
    Nf=atoi(argv[2]);
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

  d Mach=0.1;
  d kappa=0.1;
  Globalconf gc(Nf,f,"air",T0,p0,kappa);
  Globalconf air=Globalconf::airSTP(Nf,f);
  air=gc;
  Geom geom1=Geom::CylinderBlApprox(gp,L,rtube);
  // IsentropicTube t1(geom1);
  HopkinsLaminarDuct t1(geom1,gc.T0,gc.T0+10);
  var pL(gc,0);
  if(Nf>0)
    pL.set(1,1);
  tube::LeftPressure bcleft(pL);
  // tube::RightIsoTWall bcright(0,T0+10);
  // TRACE(100,"Add bc to tube...");
  t1.addBc(bcleft);
  EngineSystem sys(air);
  sys.setTimingConstraint(0,0,3,2);
  sys.setAmplitudeDof(0,0,3,1);
  // TaSystem sys(air);  
  sys.addSeg(t1); 
  sys.init();

  Solver sol(sys);

  // sol.doIter();
  // sol.solve();
    // // // vd x=t1.GetRmomes();
  // // vd err=sol.sys->Error();
  // // cout << "error:\n"<<err;

  // cout << "err:\n"<<er;
  edmat Jac=sol.sys().jac();
  // sol1.sys->show();
  cout << "Jac:\n"<< Jac<<"\n";
  // sol.solve();
  // cout << "\nDeterminant Jac:" << arma::det(Jac) << "\n";
  // sol.sys.show(true);
  return 0;
}
