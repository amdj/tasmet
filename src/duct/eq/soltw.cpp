// soltw.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////
#include "solidh.h"
#include "soltw.h"
#include "cell.h"
#include "jacrow.h"
#include "laminarduct.h"
#include "geom.h"

namespace duct {
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;
  using tasystem::var;

  SolTw::SolTw(const Cell& v,const LaminarDuct& t):Equation(v){
    TRACE(15,"SolTw::SolTw()");
    duct=&t;
    if(duct->hasSolid()){
      assert(!duct->isInsulated());
      solid=&t.getSolid();
      h=new SolidH(t.geom().shape());
    }
  }
  SolTw::~SolTw(){
    delete h;
  }
  void SolTw::init() {
    TRACE(15,"SolTw::init()");
  }
  void SolTw::show() const{
    TRACE(15,"SolTw::show()");
    cout << "------------- SolTw equation\n";
  }

  vd SolTw::error() const {		// Error in momentum equation
    TRACE(15,"SolTw::Error(), i="<<v.geti());

    if(duct->isInsulated()){
      assert(!solid);
      return v.T()()-v.Tw()();	// Set wall temp equal to fluid
    // temp. Such that dTwdx equals dTdx
    }
    else{
      assert(solid);      
      return v.vSs*h->H(v,*solid)*(v.Ts()()-v.Tw()())	\
	-duct->heatSource().Qsf(v);
    }
  }
  JacRow SolTw::jac() const {
    TRACE(15,"SolTw::jac()");
    if(duct->isInsulated()){
      assert(!solid);
      JacRow jac(dofnr,2);
      jac+=JacCol(v.T() ,eye(v));
      jac+=JacCol(v.Tw(),-eye(v));
      return jac;
    }
    else{
      assert(solid);
      JacRow jac(dofnr,6);
      dmat H=h->H(v,*solid);
      jac+=JacCol(v.Ts(),v.vSs*H);
      jac+=JacCol(v.Tw(),-v.vSs*H);
      jac+=(duct->heatSource().dQsf(v)*=-1);
      return jac;
    }
  }
  
  void SolTw::domg(vd& domg) const {
    TRACE(15,"SolTw::domg()");
  }

} // namespace duct


//////////////////////////////////////////////////////////////////////
