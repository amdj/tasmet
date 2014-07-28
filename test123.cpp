/* test.cpp */

#include "globalconf.h"
#include "tube.h"
#include "twimpedance.h"
#include "pressurebc.h"
#include "tubevertex.h"
#include "impedancebc.h"
#include "system.h"
#include "solver.h"
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
  double f=100;
  double omg=2*number_pi*f;
  double T=1/f;
  if(argc>1)
    loglevel=atoi(argv[1]);
  if(argc>2)
    Nf=atoi(argv[2]);
  cout<< "Loglevel:"<<loglevel<<"\n";
  initlog(loglevel);
  d L=0.01;
  d rtube=5e-3;
  d S=number_pi*pow(rtube,2);
  d phi=1.0;
  d PI=2*number_pi*rtube;
  d rh=S/PI;
  d p0=101325.0;
  d T0=293.15;
  d griddx=L/gp;
  d S0=S;

  d Mach=0.1;
  d kappa=0.1;
  Globalconf gc(Nf,f,"air",T0,p0,Mach,S0,griddx,0,kappa);
  Geom geom1=Geom::Cylinder(gp,L,rtube);
  Tube t1(geom1);
  t1.Init(gc);
  // var pL(gc,0);
  // if(Nf>0)
    // pL.set(1,1);
  // tube::LeftPressure bcleft(0,pL);
  // tube::RightImpedance bcright(0,415*vd(Ns,fillwith::ones));
  // tube::TwImpedance bcright(1);
  
  TAsystem sys(gc);
  sys.addseg(t1);
  sys.addseg(t1);  
  sys.connectSegs(0,1,SegCoupling::tailhead);
  // sys.addbc(bcright);
  // sys.addbc(bcleft);
  // cout << sys.getSeg(0)->GetRes() << "\n";
  Solver sol(sys);
  // sol.sys->Init();
  // sys.getSeg(0)->Error();
  // sol.sys->getSeg(0)->Jac();

  sol.sys->Init();
  // sys.show(true);
  // sol.sys->show(true);
  cout << "Jac:\n"<< sol.sys->Jac();
// sol.sys

  // Solver sol2(sol);
  
  // sol.DoIter();
    // // // vd x=t1.GetRmomes();
  // // vd err=sol.sys->Error();
  // // cout << "error:\n"<<err;



  // cout << "err:\n"<<er;
  // dmat Jac=sol.sys->Jac();
  // sol1.sys->show();
  // cout << "Jac:\n"<< Jac;

  return 0;
}
