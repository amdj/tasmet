/* test.cpp */
#include <vtypes.h>
#include "globalconf.h"
#include "tube/tube.h"
using namespace std;

int main() {
  cout <<  "Running test..." << endl;
  initlog(-1);
  us gp=11; gp=4;
  us Nf=0;
  us Ns=2*Nf+1;
  double f=100;
  double omg=2*pi*f;
  double T=1/f;
  d L=1;
  d rtube=1e-3;
  d S=pi*pow(rtube,2);
  d phi=1.0;
  d PI=2*pi*rtube;
  d rh=S/PI;
  d p0=101325.0;
  d T0=293.15;
  globalconf::Globalconf gc(Nf,f,"air");
  tube::Geom geom1(gp,L,S,phi,rh,"blapprox");
  tube::Tube t1(gc,geom1);
  t1.Init(T0,p0);
  TRACE(0,t1.vvertex[1].rho.tdata());
  
  // // TRACE(0,"vVf"<<t1.geom.vVf);
  // TRACE(0,"testoutput:"<<t1.vvertex[1].m.Error());
  // TRACE(0,t1.vvertex[    TRACE(0,t1.vvertex[1].wRr);
  TRACE(0,t1.vvertex[1].wLl);
  TRACE(0,t1.vvertex[1].wRr);
  TRACE(0,t1.vvertex[1].wLr);
  TRACE(0,t1.vvertex[1].wRl);
  TRACE(0,t1.vvertex[1].eq[0]->wLl);
  TRACE(0,t1.vvertex[1].eq[0]->wRr);
  TRACE(0,t1.vvertex[1].eq[0]->wLr);
  TRACE(0,t1.vvertex[1].eq[0]->wRl);

  TRACE(0,"test2output"<<t1.vvertex[1].m());
  // TRACE(0,"t1.gps[1].rho()");
  // TRACE(0,t1.gps[1].rho());
  // TRACE(0,t1.vvertex[1].p());
  // TRACE(0,t1.geom.vSf);

  // TRACE(0,"t1.gps[1].c.vertex.rho();"<<t1.gps[1].c.vertex.rho());
  // TRACE(0,"Tube init starting...");
  // t1.Init(293.15,101325);
  // TRACE(10,"rho:"<<t1.gps[1].rho.tdata());

  // cout <<  t1.gps[1].m.dUi();
  // TRACE(10,"dmdpip1():");
  // cout <<  t1.gps[1].m.dpip1();


  TRACE(0,"_-----------------------------------------");

  return 0;
}
