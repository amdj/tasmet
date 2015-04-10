#include "leftcell.h"

namespace tube{
  using variable::var;
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

    // vars.push_back(&mHL_);
    // Remove momentum equation from list. Put equation for Mu in
    // place of momentum eq.
    // WARN("HERE SOME EQS NEED TO BE DELETED");
    assert(eqs.find(EqType::Mom)!=eqs.end());
    // If these elements are already deleted, we do nothing
    try{
      Equation* mom=eqs.at(EqType::Mom);
      delete mom;
      TRACE(25,"MOm eq deleted");
      eqs.erase(EqType::Mom);
    }
    catch(std::out_of_range&){}
  }
  void LeftCell::show(us detailnr) const{
    cout << "------------- LeftCell ----------\n";
    Cell::show(detailnr);
  }

} // namespace tube
