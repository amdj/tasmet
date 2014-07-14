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
  int loglevel=4;
  us gp=40;
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
  Geom geom1(gp,L,S,phi,rh,"inviscid");
  Tube t1(geom1);

  var pL(gc,0);
  if(Nf>0)
    pL.set(1,1);
  tube::LeftPressure bcleft(0,pL);
  // tube::RightImpedance bcright(0,415*vd(Ns,fillwith::ones));
  tube::TwImpedance bcright(0);
  cout << "left: "<<bcleft.left <<"\n" ;
  
  TAsystem sys(gc);
  sys.addseg(t1);
  sys.addbc(bcright);
  sys.addbc(bcleft);

  Solver sol(sys);

  sol.Init();

  Solver sol1(sol);
  sol1.Init();
  vd res=sol1.sys->GetRes();
  // cout << "res:\n"<<res;
  vd err=sol1.sys->Error();
  
  // // // vd x=t1.GetRmomes();
  // // vd err=sol.sys->Error();
  // // cout << "error:\n"<<err;

  for(us i=0; i<4;i++)  
    sol1.DoIter();
  // res=sol1.sys->GetRes();
  // vd er=sol1.sys->Error();
  // cout << "res:\n"<<res;

  // cout << "err:\n"<<er;
  // dmat Jac=sol1.sys->Jac();
  // sol1.sys->show();
  // cout << "Jac:\n"<< Jac;

  return 0;
}
