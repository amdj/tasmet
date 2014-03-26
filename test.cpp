/* test.cpp */
#include <vtypes.h>
#include "globalconf.h"
#include "tube/tube.h"
#include "tube/bcvertex.h"
#include "system.h"
using namespace std;
using namespace tasystem;

int main() {
  cout <<  "Running test..." << endl;
  initlog(-2);
  us gp=11; gp=3;
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
  tasystem::Globalconf gc(Nf,f,"air");
  tube::Geom geom1(gp,L,S,phi,rh,"blapprox");
  tube::Tube t1(gc,geom1);
  // segment::Seg t1(gc);
  variable::var presLeft(t1.vop);
  presLeft.set(p0-1,0);
  TRACE(0,presLeft.tdata());
  tube::LeftPressure* bcleft=new tube::LeftPressure(t1,presLeft,T0);
  cout << bcleft->esource();
  t1.setLeftbc(bcleft);  

  t1.Init(T0,p0);
  TRACE(0,"-----------------------------------------");
  dmat j=t1.Jac();
  vd err=t1.Error();
  TRACE(0,"Error"<<endl<<err);
  TRACE(0,"Jacobian:"<<endl<<j);
  TRACE(0,"Determinant of Jacobian:"<< det(j));
  
  vd dx=-solve(j,err);
  TRACE(0,"dx:"<<endl<<dx);

  // TRACE(0,j.row(2));
  // TRACE(0,j.row(7));
  // TRACE(0,j.row(12));  
    //-------------------------------------------------    
  // TAsystem sys(gc);
  // sys.addseg(t1);
  // vd err=sys.Error();
  // vd res=sys.GetRes();
  // cout << err;
  // TRACE(0,"Result:"<<res);
  // dmat Jac=sys.Jac();
  // cout << Jac;




  TRACE(0,"_----------------------------------------");

  return 0;
}


