/* test.cpp */

#include "globalconf.h"
#include "tube.h"
#include "isentropictube.h"
#include "twimpedance.h"
#include "pressurebc.h"
#include "hopkinslaminarduct.h"
#include "tubevertex.h"
#include "impedancebc.h"
#include "isotwall.h"
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
  Globalconf gc(Nf,f,"air",T0,p0,0,kappa);
  Globalconf air=Globalconf::airSTP(Nf,f);
  air=gc;
  Geom geom1=Geom::Cylinder(gp,L,rtube);
  IsentropicTube t1(geom1);
  // HopkinsLaminarDuct t1(geom1);
  // TRACE(100,"Thisrho:");  
  TRACE(100,"Thisrho:"<< gc.gas.rho(T0,p0));
  // gc2=gc;
  // TRACE(100,"Thisrho:");   
  // TRACE(100,"Thisrho:"<< gc2.gas.rho(T0,p0));
  // TRACE(100,"Thisrho:");  

  var pL(gc,0);
  if(Nf>0)
    pL.set(1,1);
  tube::LeftPressure bcleft(0,pL);
  // tube::RightImpedance bcright(0,415*vd(Ns,fillwith::ones));
  // tube::TwImpedance bcright(0);
  // tube::RightIsoTWall bcright(0,T0+10);
  // TRACE(100,"Add bc to tube...");
  t1.addBc(bcleft);
  // TRACE(100,"Add bc to tube done...");

  // TRACE(100,"Creating empty TAsystem...");
  TAsystem sys(air);
  // TRACE(100,"Filling TAsystem with a segment...");
  // auto t2(t1);
  sys.addSeg(t1); 
  // TRACE(100,"Filling segment done...");      
  // // sys.connectSegs(0,1,SegCoupling::tailhead);
  sys.init();
  // // sys.addBc(bcright);
  // TRACE(100,"Copy TAsystem...");    
  // // TAsystem sys2(gc);
  // // sys2.addSeg(t2);
  // TRACE(100,"Copy TAsystem done...");      
  // TRACE(100,"Running init on first sys:=================================");
  // TRACE(100,"---------------------------------------------------sys init done. now starting sys2 init9");
  // sys2.init();
  // sys2.show();
  // cout << sys.getSeg(0)->GetRes() << "\n";
  // Solver sol(sys);
  // sol.sys.init();
  // sys.getSeg(0)->Error();
  // sol.sys->getSeg(0)->Jac();
  // sol.sys->init();
  // sys.show(true);
  // sol.sys->show(true);
  // cout << "Jac:\n"<< sol.sys->jac();
// sol.sys

  // Solver sol2(sol);n
  cout << "Result:\n" << sys.getRes();  
  cout << "Error:\n" << sys.error();
  // sol.doIter();
    // // // vd x=t1.GetRmomes();
  // // vd err=sol.sys->Error();
  // // cout << "error:\n"<<err;

  // cout << "err:\n"<<er;
  dmat Jac=sys.jac();
  // sol1.sys->show();
  cout << "Jac:\n"<< Jac;
  // sol.sys.show(true);
  return 0;
}

