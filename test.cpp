/* test.cpp */
#include <vtypes.h>
#include "globalconf.h"
#include "tube/tube.h"
using namespace std;

int main() {
  cout <<  "Running test..." << endl;
  initlog(-4);
  us gp=11; gp=11;
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

    // vd sinus(Ns);
    // for(us i=0;i<Ns;i++) { sinus(i)=cos(2*pi*i/Ns);}
    // vd test=init+sinus;
    globalconf::Globalconf gc(Nf,f,"air");
  tube::Geom geom1(gp,L,S,phi,rh,"blapprox");
  tube::Tube t1(gc,geom1);
  TRACE(0,"wrR 0:"<< t1.gps[9]->wRr);
  TRACE(0,"wRl 0:"<< t1.gps[9]->wRl);
  TRACE(0,"wLl 0:"<<t1.gps[1]->wLl);
  TRACE(0,"vVf"<<t1.geom.vVf);
  // TRACE(0,"t1.gps[1]->rho()");
  // TRACE(0,t1.gps[1]->rho());


  // TRACE(0,"t1.gps[1]->c.vertex->rho();"<<t1.gps[1]->c.vertex->rho());
  // TRACE(0,"Tube init starting...");
  // t1.Init(293.15,101325);
  // TRACE(10,"rho:"<<t1.gps[1]->rho.tdata());

  // cout <<  t1.gps[1]->m.dUi();
  // TRACE(10,"dmdpip1():");
  // cout <<  t1.gps[1]->m.dpip1();


  TRACE(0,"_-----------------------------------------");
  // TRACE(0,"t1.gps[1]->c.vertex->rho();"<<t1.gps[1]->c.vertex->rho());
  // t1.gps[0]->rho.set(1.2,0);
  // cout << t1.gps[1]->c.Error();
  // cout << t1.gps[1]->m.Error();
  // cout << t1.gps[0]->s.Error();
  // cout << t1.gps[1]->e.Error();
  // TRACE(10,"dEnergy/dpim1:"<< t1.gps[1]->e.dpim1());
  // TRACE(10,"Momentum:dpip1"<< t1.gps[1]->e.dpip1());


  // TRACE(10,"Error at gp 1:"<< t1.gps[1]->Error());
  // // TRACE(0,"Gridpoint 1 Jacobian:"<<t1.gps[1]->Jacobian());
  // // t1.gps[1]->Jacobian();
  // vd res1=t1.gps[1]->Get();
  // // res1(1)=3333;
  // t1.gps[1]->Set(res1);

  // vd totres=t1.Get();
  // totres(1)=555;

  // t1.Set(totres);
  // TRACE(0,"U(0):"<<t1.gps[0]->U());
  // totres=t1.Get();
  // TRACE(10,"Total result:"<<totres);
  // TRACE(10,"Total error:"<<t1.Error());
  return 0;
}
