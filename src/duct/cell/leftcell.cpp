#include "leftcell.h"

namespace duct{
  using tasystem::var;
  using tasystem::JacRow;
  using tasystem::JacCol;


  
  LeftCell::LeftCell(us i,const Duct& t):
    BcCell(i,t)
  {
    TRACE(15,"LeftCell::LeftCell()");
  }
  void LeftCell::init(const Cell* left,const Cell* right){
    TRACE(10,"LeftCell::init()");
    assert(!left);
    assert(right);
    BcCell::init(left,right);
  }
  void LeftCell::show(us detailnr) const{
    cout << "------------- LeftCell ----------\n";
    Cell::show(detailnr);
  }
  void LeftCell::updateNf(){
    TRACE(15,"LeftCell::updateNf()");
    BcCell::updateNf();
  }
} // namespace duct
