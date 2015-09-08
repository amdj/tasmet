// solidh.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////
#include "var.h"
#include "solidh.h"
#include "cell.h"
#include "solid.h"
#include "exception.h"

namespace tube {
  using tasystem::var;
  using solids::Solid;

  class SolidHVert:public SolidH{

  public:
    SolidHVert():SolidH(){
    }
    dmat H(const Cell& v,const Solid& s) const {
      var H(*v.gc);

      // Since Ss=b*2*rhs
      // And Sf=b*2*rh
      // Divide the two to find:
      // Ss/Sf=rhs/rh
      //
      // So: rhs=rh*Ss/Sf
      d rhs=v.vrh*v.vSs/v.vSf;	// Hydraulic radius of solid
      d kappas=s.kappa(v.Ts()(0));
      H.setadata(0,3*kappas/pow(rhs,2));

      return H.freqMultiplyMat();
    }
  };
  SolidH::SolidH(const string& shape){
    if(shape.compare("vert")==0)
      h=new SolidHVert();
    else
      throw MyError("Modeling the solid for shapes other than"
		    " VertPlates is not yet implemented.");
  }

} // namespace tube
//////////////////////////////////////////////////////////////////////
