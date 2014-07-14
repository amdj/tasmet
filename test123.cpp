/* test.cpp */

#include "globalconf.h"
#include "tube/tube.h"
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
  d kappa=0;
  Globalconf gc(Nf,f,"air",T0,p0,Mach,S0,griddx,0,kappa);
  variable::var U(gc);
  // U.set(1e-2,1);

  U.set(1e0,3);    
  // U.set(1e0,2);
  TRACE(10,"U:"<<gc.fDFT*(U.tdata()));  
  TRACE(10,endl<<gc.iDFT);
  TRACE(10,"Usq::"<<(U*U)());

  // TRACE(10,gc.DDTfd);

  Geom geom1(gp,L,S,phi,rh,"inviscid");
  Tube t1(geom1);
  // TRACE(10,"Dl:"<<endl<<t1.vvertex[1]->eq[0]->D_l());
  // TRACE(10,"-Dr:"<<endl<<-t1.vvertex[1]->eq[0]->D_r());  

  // vd Z=(415/S)*ones<vd>(gc.Ns);
  // tube::RightImpedance bcright(geom1,Z);
  var pL(gc,0);
  pL.set(1,1);
  tube::LeftPressure bcleft(0,pL);
  tube::RightImpedance bcright(0,415*vd(Ns,fillwith::ones));
  
  TAsystem sys(gc);
  sys.addseg(t1);
  sys.addbc(bcright);
  sys.addbc(bcleft);
  // // // // TRACE(0,bcright->Z);  
  Solver sol(sys);
  sol.Init();

  Solver sol1(sol);
  // cout << "Ndofs for seg before adding :"<< t1.getNdofs()<< "\n";
  // cout << "Ndofs for seg in system :"<< sys[0]->getNdofs()<< "\n";
  // cout << "Ndofs for seg in solver:"<< (*sol.sys)[0]->getNdofs()<< "\n";

  // cout << "Ncells for seg before adding :"<< t1.getNcells()<< "\n";
  // cout << "Ncells for seg in system :"<< sys[0]->getNcells()<< "\n";
  // cout << "Ncells for seg :"<< (*sol.sys)[0]->getNcells()<< "\n";

  // sol1.Init();

  vd res=sol1.sys->GetRes();
  cout << "res:\n"<<res;
  // vd err=sol1.sys->Error();
  
  // // vd x=t1.GetRmomes();
  // vd err=sol.sys->Error();
  // cout << "error:\n"<<err;

  // for(us i=0; i<4;i++)  
    // sol1.DoIter();
  res=sol1.sys->GetRes();
  vd er=sol1.sys->Error();
  cout << "res:\n"<<res;

  cout << "err:\n"<<er;
  dmat Jac=sol1.sys->Jac();
  sol.sys->show();
  cout << "Jac:\n"<< Jac;

  return 0;
}
