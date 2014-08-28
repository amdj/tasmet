#include "solver.h"
#include "var.h"
#include "pressurebc.h"
#include "impedancebc.h"

#include "seg.h"
#include "tube.h"
#include "geom.h"
using namespace tasystem;
typedef std::shared_ptr<Solver> Solverptr;
typedef std::shared_ptr<tube::Geom> Geomptr;
typedef std::shared_ptr<tube::Tube> Tubeptr;



Solverptr TwoTubes(int loglevel=10)
{
  us Nf=5;

  using segment::vertexptr;
  using math_common::esdmat;
  using namespace segment;
  using tube::Tube;

  inittrace(loglevel);

  us Ns=2*Nf+1;
  d freq=100;
  d L=6.0/2; 			// Two tubes of half the length
  d S=0.1;
  d R=sqrt(S/number_pi);
  d rh=R/2;
  d PI=S/rh;
  us gp=121/2;
   // gp=4;
  Globalconf gc(Nf,freq);
  segment::Geom g1(gp,L,S,1.0,S/PI,"circ");
  Tube t1(g1);
  Tube t2(g1);
  d p1=1;
  variable::var presLeft(gc);
  if (Nf>0)
    presLeft.set(p1,1);	// One oscillation
  vertexptr bcleft1(new tube::LeftPressure(t1,presLeft));
  t1.setLeftbc(bcleft1);
  // vertexptr bcleft2(new tube::LeftPressure(t2,presLeft));
  // t2.setLeftbc(bcleft2);

  vd Z=(415/S)*ones<vd>(Ns);
  // vertexptr bcright1(new tube::RightImpedance(t1,Z));
  vertexptr bcright2(new tube::RightImpedance(t2,Z));
  // t1.setRightbc(bcright1);
  t2.setRightbc(bcright2);
  coupleSegs(t1,t2,tailhead);
  TaSystemptr sys(new tasystem::TaSystem(gc));
  sys->addseg(t1);
  sys->addseg(t2);
  Solverptr sol(new tasystem::Solver(*sys));

  return sol;
}













