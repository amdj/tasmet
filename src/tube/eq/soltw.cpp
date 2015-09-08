// soltw.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////
#define TRACERPLUS 60
#include "solidh.h"
#include "soltw.h"
#include "cell.h"
#include "jacrow.h"
#include "laminarduct.h"
#include "geom.h"

namespace tube {
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;
  using tasystem::var;

  SolTw::SolTw(const Cell& v,const LaminarDuct& t):Equation(v){
    TRACE(15,"SolTw::SolTw()");
    solid=&t.getSolid();
    tube=&t;
    h=new SolidH(t.geom().shape());
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
    assert(solid);
    return v.vSs*h->H(v,*solid)*(v.Ts()()-v.Tw()())\
      -tube->heatSource().Qsf(v);
  }
  JacRow SolTw::jac() const {
    TRACE(15,"SolTw::jac()");
    assert(solid);
    JacRow jac(dofnr,6);
    jac+=JacCol(v.Ts(),v.vSs*h->H(v,*solid));
    jac+=JacCol(v.Tw(),-v.vSs*h->H(v,*solid));
    jac+=(tube->heatSource().dQsf(v)*=-1);
    return jac;
  }
  
  void SolTw::domg(vd& domg) const {
    TRACE(15,"SolTw::domg()");
  }

} // namespace tube


//////////////////////////////////////////////////////////////////////
