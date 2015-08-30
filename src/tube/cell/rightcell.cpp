#include "rightcell.h"

namespace tube{
  using tasystem::var;
  using tasystem::JacRow;
  using tasystem::JacCol;


  // These functions should stay internal to this unit

 RightCell::RightCell(us i,const Tube& t):
    BcCell(i,t)
  {
    const tasystem::Globalconf& gc=*(this->gc);
    // Initialize right wall variables
  }
  void RightCell::init(const Cell* left,const Cell* right){
    TRACE(15,"RightCell::init()");
    assert(!right);             // Otherwise, this is not the
                                // rightmost!
    assert(left);
    BcCell::init(left,right);

    mr_=var(*gc);
    vars.push_back(&mr_);
  }
  void RightCell::show(us detailnr) const{
    cout << "------------- RightCell ---------\n";
    Cell::show(detailnr);
  }

} // namespace tube


