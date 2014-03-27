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
  initlog(15);
  us gp=100;
  us Nf=5;
  us Ns=2*Nf+1;
  double f=10;
  double omg=2*pi*f;
  double T=1/f;
  d L=0.01;
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
  presLeft.set(p0,0);
  if (Nf>0)
    presLeft.set(1,1);	// One oscillation
  TRACE(10,presLeft.tdata());
  tube::LeftPressure* bcleft=new tube::LeftPressure(t1,presLeft,T0);
  // cout << bcleft->esource();
  t1.setLeftbc(bcleft);  

  t1.Init(T0,p0);
  TRACE(10,"-----------------------------------------");

  // dmat j=t1.Jac();
  
  
  // TRACE(10,"Jacobian:"<<endl<<j);
  // TRACE(10,"Determinant of Jacobian:"<< det(j));
  



  // TRACE(10,j.row(2));
  // TRACE(10,j.row(7));
  // TRACE(10,j.row(12));  
    //-------------------------------------------------    
  TAsystem sys(gc);
  sys.addseg(t1);
  vd err=sys.Error();
  dmat j=sys.Jac();


  vd dx=-solve(j,err);
  vd xold=sys.GetRes();
  vd xnew=xold+dx;
  sys.SetRes(xnew);
  TRACE(10,"err:"<<endl<<err);
  
  TRACE(10,"xold:"<<endl<<xold);
  TRACE(10,"dx:"<<endl<<dx);  
  TRACE(10,"xnew:"<<endl<<xnew);

  TRACE(10,"-----------------------------------------");

  return 0;
}


