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
  }
  void LeftCell::show(us detailnr) const{
    cout << "------------- LeftCell ----------\n";
    Cell::show(detailnr);
  }

} // namespace tube
