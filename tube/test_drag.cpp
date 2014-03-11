#include "drag.h"
#include "../globalconf.h"
#include "tube.h"
int main(){
  cout <<  "Running test..." << endl;
  initlog(-1);

  us gp=11; gp=3;
  us Nf=3;
  us Ns=2*Nf+1;
  double f=100;
  double omg=2*pi*f;
  double T=1/f;
  // vd sinus(Ns);
  // for(us i=0;i<Ns;i++) { sinus(i)=cos(2*pi*i/Ns);}
  // vd test=init+snus;
  d r=2e-3;
  d Sf=pi*pow(r,2);
  d L=1.0;
  d phi=1.0;
  globalconf::Globalconf gc(Nf,f,"air");
  tube::Geom geom1(gp,L,Sf,phi,Sf/(2*pi*r),"vert");
  tube::Tube t1(gc,geom1);
  TRACE(0,"Geom.shape:"<< geom1.shape);
  t1.Init(293.15,101325);
  tube::LaminarDragResistance l(t1);
  vc rh_o_dnu(2,fillwith::zeros);

  rh_o_dnu(1)=1;
  t1.gps[0]->U.set(1.0,0);
  t1.gps[0]->U.set(1.0,1);
  t1.gps[0]->U.set(1.0,3);
  t1.gps[0]->U.set(1.0,5);

  TRACE(0,"Temp:"<<  t1.gps[0]->T(0));

  vd Drag=l(0);
  dmat dUi=l.dUi(0);
  TRACE(0,"dUi:"<<dUi);
  TRACE(0,"Drag:"<<Drag);
  TRACE(0,"U:"<<t1.gps[0]->U());
  
  return 0;
}
