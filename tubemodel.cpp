#include "solver.h"
#include "var/var.h"
#include "tube/bcvertex.h"
int main(int argc, char *argv[])
{
  SPOILNAMESPACE
    using segment::vertexptr;
    using math_common::esdmat;
    initlog(10);
  us Nf=0;
  us Ns=2*Nf+1;
  d freq=100;
  d L=6;
  d S=0.1;
  d R=sqrt(S/number_pi);
  d rh=R/2;
  d PI=S/rh;
  // us gp=120;
  us gp=4;
  tasystem::Globalconf gc(Nf,freq);
  tube::Geom g1(gp,L,S,1.0,S/PI,"circ");
  tube::Tube t1(gc,g1);
  variable::var presLeft(t1.gc);
  if (Nf>0)
    presLeft.set(1.0,1);	// One oscillation

  vertexptr bcleft(new tube::LeftPressure(t1,presLeft));
  t1.setLeftbc(bcleft);

  vd Z=(415/S)*ones<vd>(Ns);
  // tube::RightImpedance* bcright=new tube::RightImpedance(t1,Z);
  // t1.setRightbc(bcright);

  tasystem::TAsystem sys(gc);
  sys.addseg(t1);
  tasystem::Solver sol(sys);
  sol.DoIter();
  // esdmat J=sys.Jac();
  // dmat J=dmat(math_common::EigenToArma(sys.Jac()));
  // TRACE(10,"J:\n"<<J);
  // vd er=sys.Error();
  // TRACE(10,"Err:\n"<< er);
  // vd x=sys.GetRes();
  // vd dx=-arma::solve(J,er);
  return 0;
}













