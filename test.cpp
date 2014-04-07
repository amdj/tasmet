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
  initlog(0);
  us gp=3;
  us Nf=1;
  us Ns=2*Nf+1;
  double f=10;
  double omg=2*number_pi*f;
  double T=1/f;
  d L=0.01;
  d rtube=1e-3;
  d S=number_pi*pow(rtube,2);
  d phi=1.0;
  d PI=2*number_pi*rtube;
  d rh=S/PI;
  d p0=101325.0;
  d T0=293.15;
  tasystem::Globalconf gc(Nf,f,"air",T0,p0);
  tube::Geom geom1(gp,L,S,phi,rh,"blapprox");
  tube::Tube t1(gc,geom1);
  // segment::Seg t1(gc);
  variable::var presLeft(t1.vop);
  if (Nf>0)
    presLeft.set(1.0,1);	// One oscillation
  // TRACE(10,presLeft.tdata());
  tube::LeftPressure* bcleft=new tube::LeftPressure(t1,presLeft);
  vd Z=(415/S)*ones<vd>(Ns);
  tube::RightImpedance* bcright=new tube::RightImpedance(t1,Z);

  TRACE(0,"bcleftTL:"<<bcleft->TL());
  TRACE(0,"bcleftpL:"<<bcleft->pL());  

  t1.setRightbc(bcright);
  t1.setLeftbc(bcleft);

  // TRACE(0,bcright->Z);  
  vd x=t1.GetRes();
  vd errold=t1.Error();
  // TRACE(0,"rho(0)"<<t1.gas.Rs()*t1.vvertex[0]->rho.tdata()%t1.vvertex[0]->T.tdata());
  
  // TRACE(10,"-----------------------------------------");
  dmat jac=dmat(t1.Jac());
  // t1.DoIter();
  // t1.DoIter();  

  // vd errnew=t1.Error();
  TRACE(0,"Jac:"<<endl<<jac);
  // dmat Jace1=jac.rows(2*Ns,2*Ns+Ns-1);
  // TRACE(10,endl<<Jace1.cols(0,5*Ns-1));
  // TRACE(10,endl<<Jace1.cols(5*Ns,10*Ns-1));
  TRACE(10,"Previous error:"<< endl<< errold);
  // TRACE(10,"New error:"<<endl<<errnew);
  TRACE(10,"xold:"<<endl<<x);
  // x=t1.GetRes();
  TRACE(10,"xnew:"<<endl<<x);
  // dmat Jace2=jac.rows(7*Ns,7*Ns+Ns-1);
  // TRACE(10,endl<<Jace2.cols(0,5*Ns-1));
  // TRACE(10,endl<<Jace2.cols(5*Ns,10*Ns-1));
  // TRACE(0,"Error:"<<endl<<err);
  // // TRACE(0,bcleft->TL());
  // TRACE(10,"-----------------------------------------");
  // TRACE(0,"Momentum::dUi() van gp 2:"<<endl<< t1.vvertex[1]->eq[1]->Jac());

  return 0;
}




