/* test.cpp */

#include "globalconf.h"
#include "tube/tube.h"
#include "tube/bcvertex.h"

using namespace std;
using namespace tasystem;

int main() {
  cout <<  "Running test..." << endl;
  initlog(2);
  us gp=4;
  us Nf=2;
  us Ns=2*Nf+1;
  double f=100;
  double omg=2*number_pi*f;
  double T=1/f;
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
  d kappa=0;
  tasystem::Globalconf gc(Nf,f,"air",T0,p0,Mach,S0,griddx,0,kappa);

  variable::var U(gc);
  // U.set(1e-2,1);

  U.set(1e0,3);    
  // U.set(1e0,2);
  TRACE(10,"U:"<<gc.fDFT*(U.tdata()));  
  TRACE(10,endl<<gc.iDFT);
  TRACE(10,"Usq::"<<(U*U)());

  TRACE(10,gc.DDTfd);
  // tube::Geom geom1(gp,L,S,phi,rh,"inviscid");
  // tube::Tube t1(gc,geom1);

  // variable::var presLeft(t1.gc);
  // if (Nf>0)
  //   presLeft.set(1.0,1);	// One oscillation
  // TRACE(10,presLeft.tdata());
  // tube::LeftPressure* bcleft=new tube::LeftPressure(t1,presLeft);
  // vd Z=(415/S)*ones<vd>(Ns);


  // t1.setLeftbc(bcleft);

  // TRACE(10,"Dl:"<<endl<<t1.vvertex[1]->eq[0]->D_l());
  // TRACE(10,"-Dr:"<<endl<<-t1.vvertex[1]->eq[0]->D_r());  
  // // tube::RightImpedance* bcright=new tube::RightImpedance(t1,Z);
  // // t1.setRightbc(bcright);
  // // // TRACE(0,bcright->Z);  

  // vd x=t1.GetRmomes();
  // vd er=t1.Error();
  // TRACE(10,"-----------------------------------------");
  // dmat jac=dmat(t1.Jac());
  // for(us h=0;h<2;h++)
    // t1.DoIter();
  // vd dx=-solve(jac,er);
  
  // vd errnew=t1.Error();
  // jac=dmat(t1.Jac());
  // dx+=-solve(jac,er);  
  // TRACE(10,"Previous error:"<< endl<< errold);
  // TRACE(10,"New error:"<<endl<<errnew);
  // TRACE(10,"xold:"<<endl<<x);
  // TRACE(0,"Jac:"<<endl<<jac);


  // TRACE(10,"xnew:"<<endl<<x+dx);

  // TRACE(10,t1.vvertex[1]->eq[1]->dxp);
  // TRACE(10,t1.vvertex[1]->eq[1]->dxm);
  // TRACE(10,gc.dx);
  // TRACE(10,"-----------------------------------------");
  // TRACE(8,t1.vvertex[0]->eq[0]->D_l());
  // TRACE(0,t1.vvertex[0]->left); 
  // TRACE(0,t1.vvertex[t1.Ncells-1]->right); 
  return 0;
}

