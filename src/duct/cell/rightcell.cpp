#include "rightcell.h"

namespace duct{
  using tasystem::var;
  using tasystem::JacRow;
  using tasystem::JacCol;


  // These functions should stay internal to this unit

 RightCell::RightCell(us i,const Duct& t):
    BcCell(i,t)
  {
    mr_=var(*gc);
  }
  void RightCell::init(const Cell* left,const Cell* right){
    TRACE(15,"RightCell::init()");
    assert(!right);             // Otherwise, this is not the
                                // rightmost!
    assert(left);
    BcCell::init(left,right);
    vars.push_back(&mr_);
  }
  void RightCell::show(us detailnr) const{
    cout << "------------- RightCell ---------\n";
    Cell::show(detailnr);
  }
  void RightCell::updateNf() {
    TRACE(15,"RightCell::updateNf()");
    BcCell::updateNf();
    mr_.updateNf();
  }
} // namespace duct


