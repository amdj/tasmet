#include "solver.h"
#include "var/var.h"
#include "bcvertex.h"
#include "seg.h"
SPOILNAMESPACE
int main(int argc, char *argv[])
{
  us Nf=5;
  us loglevel=0;
  if(argc>1)
    loglevel=atoi(argv[1]);
  if(argc>2)
    Nf=atoi(argv[2]);
  cout<< "Loglevel:"<<loglevel<<"\n";

  using segment::vertexptr;
  using math_common::esdmat;
  using namespace segment;
  initlog(loglevel);

  us Ns=2*Nf+1;
  d freq=100;
  d L=6.0/2; 			// Two tubes of half the length
  d S=0.1;
  d R=sqrt(S/number_pi);
  d rh=R/2;
  d PI=S/rh;
  us gp=120/2;
   // gp=4;
  tasystem::Globalconf gc(Nf,freq);
  tube::Geom g1(gp,L,S,1.0,S/PI,"circ");
  tube::Tube t1(gc,g1);
  tube::Tube t2(gc,g1);
  
  coupleSegs(t1,t2,tailhead);
  variable::var presLeft(t1.gc);
  if (Nf>0)
    presLeft.set(100.0,1);	// One oscillation

  vertexptr bcleft1(new tube::LeftPressure(t1,presLeft));
  t1.setLeftbc(bcleft1);
  vertexptr bcleft2(new tube::LeftPressure(t2,presLeft));
  t2.setLeftbc(bcleft2);

  vd Z=(415/S)*ones<vd>(Ns);
  vertexptr bcright1(new tube::RightImpedance(t1,Z));
  vertexptr bcright2(new tube::RightImpedance(t2,Z));
  t1.setRightbc(bcright1);
  t2.setRightbc(bcright2);
  tasystem::TAsystem sys(gc);
  sys.addseg(t1);
  sys.addseg(t2);
  tasystem::Solver sol(sys);
  TRACE(15,"Error:"<< arma::norm(sys.Error(),2));
  sol.DoIter();
  TRACE(15,"Error:"<< arma::norm(sys.Error(),2));
  sol.DoIter();
  TRACE(15,"Error:"<< arma::norm(sys.Error(),2));

  return 0;
}













