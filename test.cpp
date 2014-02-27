/* test.cpp */
#include "common/vtypes.h"
#include "globalconf.h"
#include "tube.h"
using namespace std;

int main() {
  cout <<  "Running test..." << endl;
  initlog(-1);
  us gp=10;
  us Nf=1;
  us Ns=2*Nf+1;
  double f=100;
  double omg=2*pi*f;
  double T=1/f;
  double init=10.0;
  // vd sinus(Ns);
  // for(us i=0;i<Ns;i++) { sinus(i)=cos(2*pi*i/Ns);}
  // vd test=init+sinus;
  globalconf::Globalconf gc(Nf,f,"air");
  tube::Geom geom1(5,1.0,0.1,1.0,sqrt(0.1/pi)/2,"circ");
  tube::Tube t1(gc,geom1);
  t1.Init(293.15,101325);
  cout <<  t1.gps[1].lc();
  d k=omg/343;
  vd ugok(gp),pgok(gp);

  return 0;
}
