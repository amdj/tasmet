#include "leftcell.h"

namespace tube{
  using tasystem::var;
  using tasystem::JacRow;
  using tasystem::JacCol;


  
  LeftCell::LeftCell(us i,const Tube& t):
    BcCell(i,t)
  {
    TRACE(15,"LeftCell::LeftCell()");
  }
  void LeftCell::init(const Cell* left,const Cell* right){
    TRACE(10,"LeftCell::init()");
    BcCell::init(left,right);
    assert(!left);
    assert(right);

    // Remove momentum equation from list. Put equation for Mu in
    // place of momentum eq.
    assert(eqs.find(EqType::Mom)!=eqs.end());
    Equation* mom=eqs.at(EqType::Mom);
    delete mom;
    TRACE(25,"MOm eq deleted");
    eqs.erase(EqType::Mom);

  }
  void LeftCell::show(us detailnr) const{
    cout << "------------- LeftCell ----------\n";
    Cell::show(detailnr);
  }

} // namespace tube
