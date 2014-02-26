/* test_lintube.cpp */
#include "common/vtypes.h"
#include "globalconf.h"
using namespace std;

int main() {
  cout <<  "Running test..." << endl;
  us gp=10;
  us Nf=2;
  us Ns=2*Nf+1;
  double f=100;
  double omg=2*pi*f;
  double T=1/f;
  double init=10.0;
  // vd sinus(Ns);
  // for(us i=0;i<Ns;i++) { sinus(i)=cos(2*pi*i/Ns);}
  // vd test=init+sinus;
  

  d k=omg/343;
  vd ugok(gp),pgok(gp);
  //  vd x=lt.getx();
  //  d L=1.0;
  //  d up=1.0;


  return 0;
}
