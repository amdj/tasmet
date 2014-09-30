/* enginetest.cpp */

#include "tube.h"
#include "isentropictube.h"
#include "hopkinslaminarduct.h"
#include "enginesystem.h"
#include "solver.h"
#include "bc.h"

using namespace std;
using namespace segment; 
using namespace tasystem;
using namespace tube;

us gp=4;
us Nf=1;
us Ns=2*Nf+1;
double f=100;
double omg=2*number_pi*f;
double T=1/f;
d p0=101325.0;
d T0=293.15;

d Tr=T0;

Solver drivenres(int argc,char* argv[]){

  if(argc>2)
    Nf=atoi(argv[2]);
  if(argc>3)
    gp=atoi(argv[3]);
  d L=0.01;

  d S=1;
  
  d rtube=sqrt(S/number_pi);

  d phi=1.0;
  d PI=2*number_pi*rtube;
  d rh=S/PI;
  d griddx=L/gp;
  d S0=S;
  d kappa=0.1;
  Globalconf gc=Globalconf::airSTP(Nf,f);
  Geom geom1=Geom::CylinderBlApprox(gp,L,rtube);
  HopkinsLaminarDuct t1(geom1,T0,Tr);
  // IsentropicTube t1(geom1);
  var pL(gc,0);
  // pL.set(0,3.14);
  if(Nf>0)
    pL.set(1,1);
  tube::LeftPressure bcleft(pL);
  // tube::RightAdiabaticWall raw;
  // tube::TwImpedance raw;  
  tube::RightIsoTWall bcright(T0+10);
  // TRACE(100,"Add bc to tube...");
  t1.addBc(bcleft);
  t1.addBc(bcright);  
  TaSystem sys(gc);  
  sys.addSeg(t1); 
  sys.init();
  return Solver(sys);

}



int main(int argc,char* argv[]) {
  cout <<  "Running enginetest..." << endl;

  Solver drivensol=drivenres(argc,argv);
  drivensol.solve();
  evd res=drivensol.sys().getRes();
  int loglevel=25;
  if(argc>1)
    loglevel=atoi(argv[1]);
  cout<< "Loglevel:"<<loglevel<<"\n";
  inittrace(loglevel);

  // cout << "Result:\n"<< res << "\n";
  TaSystem& sys=drivensol.sys();
  // Tube& t=asTube(*sys.getSeg(0));
  Globalconf gc=sys.gc;
  gc.setDriven(false);


  Geom geom1=drivensol.sys().getSeg(0)->geom;
  // IsentropicTube t1(geom1);
  HopkinsLaminarDuct t1(geom1,gc.T0);
  t1.addBc(LeftAdiabaticWall());
  // t1.addBc(LeftEnginePressure(1));
  t1.addBc(RightIsoTWall(gc.T0+10));
  EngineSystem esys(gc);
  // TaSystem esys(gc);
  // esys.setTimingConstraint(0,0,2,2);
  // esys.setAmplitudeDof(0,0,2,1);
  
  esys.addSeg(t1);
  Solver newsol(esys);
  newsol.sys().init();
  newsol.sys().setRes(sys);
  // asTube(*newsol.sys().getSeg(0)).setRes(t);
  // esys.setRes(res);
  cout << "error:\n" << newsol.sys().error() << "\n";
  // vd domg=dynamic_cast<EngineSystem&>(newsol.sys()).domg();
  // cout << "domg: "<< domg << "\n";

  // cout << "dmtotdx:" << dmtotdx << "\n";
  // drivensol.sys().showJac();
  newsol.sys().showJac();
  // newsol.doIter();
  // newsol.solve();
  return 0;
}

