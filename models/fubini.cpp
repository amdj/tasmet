#include "models.h"


#include "gas.h"
#include "isentropictube.h"
#include "hopkinslaminarduct.h"
#include "bc.h"

using namespace segment;
using namespace tasystem;
using namespace tube;
using namespace variable;
using namespace gases;
Solver* Fubini(us gp,us Nf,d freq,d L,vd p1,int loglevel,d kappa,int options)
{
  inittrace(loglevel);
  d dx=L/gp;
  d S=1;
  d T0=293.15;
  d p0=101325;
  Gas g("air");
  d rho0=g.rho(T0,p0);
  d c0=g.cm(T0);
  d z0=rho0*c0;
  d R=sqrt(S/number_pi);

  d PI=S/(2*number_pi*R);
  d Mass=0;
  us Ns=2*Nf+1;
  cout << "Kappa: " << kappa << "\n";
  Globalconf gc=Globalconf::airSTP(Nf,freq);
  gc.show();
  Geom geom1=Geom::CylinderBlApprox(gp,L,R);

  // TRACE(30,"p1:"<<p1);
  var pL(gc);
  for(us i=0;i<Ns;i++)
    pL.set(i,p1(i));

  // vd Zv=(z0/S)*vd(2*Nf+1,fillwith::ones);

  LeftPressure bleft(pL);
  TwImpedance bright;
  // RightImpedance bright(0,Zv);
  // RightIsoTWall bright(0,T0);
  TaSystem sys(gc);

  if(options & ISENTROPIC)  {
    IsentropicTube t1(geom1);
    t1.addBc(bleft);
    t1.addBc(bright);
    sys.addSeg(t1);
  }
  else{
    HopkinsLaminarDuct t1(geom1,gc.T0);
    t1.addBc(bleft);
    t1.addBc(bright);
    sys.addSeg(t1);
  }

  Solver* Sol=new Solver(sys);
  Sol->sys().show(false);
  return Sol;  
}



