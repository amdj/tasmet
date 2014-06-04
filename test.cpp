/* test.cpp */

#include "globalconf.h"
#include "tube/tube.h"
#include "pressurebc.h"
#include "tubevertex.h"
#include "impedancebc.h"
#include "system.h"
#include "solver.h"
using namespace std;
using namespace tasystem;

int main(int argc,char* argv[]) {
  cout <<  "Running test..." << endl;
  us loglevel=10;
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
  d kappa=0.25;kappa=0;
  tasystem::Globalconf gc(Nf,f,"air",T0,p0,Mach,S0,griddx,0,kappa);

  variable::var U(gc);
  // U.set(1e-2,1);

  U.set(1e0,3);    
  // U.set(1e0,2);
  TRACE(10,"U:"<<gc.fDFT*(U.tdata()));  
  TRACE(10,endl<<gc.iDFT);
  TRACE(10,"Usq::"<<(U*U)());

  // TRACE(10,gc.DDTfd);
  tube::Geom geom1(gp,L,S,phi,rh,"inviscid");
  tube::Tube t1(gc,geom1);

  variable::var presLeft(t1.gc);
  if (Nf>0)
    presLeft.set(1.0,1);	// One oscillation
  // TRACE(10,presLeft.tdata());
  tube::LeftPressure* bcleft=new tube::LeftPressure(t1,presLeft);
  vd Z=(415/S)*ones<vd>(Ns);


  t1.setLeftbc(bcleft);
  // TRACE(10,"Dl:"<<endl<<t1.vvertex[1]->eq[0]->D_l());
  // TRACE(10,"-Dr:"<<endl<<-t1.vvertex[1]->eq[0]->D_r());  
  tube::RightImpedance* bcright=new tube::RightImpedance(t1,Z);
  // t1.setRightbc(bcright);
  tasystem::TAsystem sys(gc);
  sys.addseg(t1);

  // // // TRACE(0,bcright->Z);  

  // vd x=t1.GetRmomes();
  // vd er=t1.Error();
  // TRACE(10,"-----------------------------------------");
  dmat jac=sys.Jac();
  dmat fjac=t1.vvertex[0]->Jac();
  dmat ljac=t1.vvertex[2]->Jac();

  // TRACE(20,"wRl:"<< ((tube::TubeVertex*) t1.vvertex[1].get())->wRl);
  // for(us h=0;h<2;h++)
    // t1.DoIter();
  // vd dx=-solve(jac,er);
  
  vd err=sys.Error();
  tasystem::Solver sol(sys);
  sol.DoIter();
  // TRACE(10,"detJ:\n"<<arma::det(jac));
  
  // vd dx=-arma::solve(jac,err);  
  // TRACE(10,"Previous error:"<< endl<< err);
  // TRACE(10,"New error:"<<endl<<errnew);
  // TRACE(10,"xold:"<<endl<<x);
  // TRACE(5,"Jac:"<<endl<<jac);
  // TRACE(5,"First Jac:"<<endl<<fjac);
  // TRACE(5,"Last Jac:"<<endl<<ljac);

  
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

