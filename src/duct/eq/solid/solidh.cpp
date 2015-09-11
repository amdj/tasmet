// solidh.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
// Implementation of computation of heat transfer coefficient of the
// solid for different geometries. Only implemented for parallel
// plates. For parallel plates, the heat transfer coefficient is
// 
//        kappa_s * i * s_s * f_s
// H = ---------------------------
//        r_h^2 * (1 - f_s)
//
//  For which the time-averaged limit is
// 
//            3 * kappa_s
// H0 = ----------------
//              r_h^2
// 
//////////////////////////////////////////////////////////////////////
#include "var.h"
#include "solidh.h"
#include "rottfuncs.h"
#include "cell.h"
#include "solid.h"
#include "exception.h"

namespace duct {
  using tasystem::var;
  using solids::Solid;

  class SolidHVert:public SolidH{

  public:
    SolidHVert():SolidH(){
    }
    dmat H(const Cell& v,const Solid& s) const {
      TRACE(15,"SolidHVert::H(Cell,solid)");
      us Nf=v.gc->Nf();
      vc H(Nf+1);

      // Since Ss=b*2*rhs
      // And Sf=b*2*rh
      // Divide the two to find:
      // Ss/Sf=rhs/rh
      //
      // So: rhs=rh*Ss/Sf
      d rhs=v.vrh*v.vSs/v.vSf;	// Hydraulic radius of solid
      d kappas=s.kappa(v.Ts()(0));

      H(0)=3*kappas/pow(rhs,2);
      
      if(Nf>0) {
	d rhoc=s.rho(v.Ts()(0))*s.c(v.Ts()(0));
	vd omgvec=v.gc->getomg()*linspace(1,Nf,Nf);

	// Solid thermal wave number
	vc  s_s=(1.0+0.0*I)*rhs*sqrt(rhoc*omgvec/kappas);
	vc  rh_ov_deltas=s_s/rottfuncs::sq2;
	// Thermal rott functi<on in solid
	vc fs=rottfuncs::f_vert(rh_ov_deltas);
	H.subvec(1,Nf)=kappas*I*pow(s_s,2)%fs/(pow(rhs,2)*(1.0-fs));
      }
    return var(*v.gc,H).freqMultiplyMat();
    }
  };
  SolidH::SolidH(const string& shape){
    TRACE(15,"SolidH::SolidH()");
    if(shape.compare("vert")==0)
      h=new SolidHVert();
    else
      throw MyError("Modeling the solid for shapes other than"
		    " VertPlates is not yet implemented.");
  }

} // namespace duct
//////////////////////////////////////////////////////////////////////
