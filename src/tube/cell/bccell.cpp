#include "bccell.h"
#include "var.h"
#include "jacrow.h"
#include "weightfactors.h"


namespace tube{
  using tasystem::JacRow;
  using tasystem::JacCol;
  using variable::var;
  BcCell::BcCell(us i,const Tube& t):
    Cell(i,t)
  {}
  void BcCell::init(const Cell* left,const Cell* right){
    TRACE(10,"BcCell::init(const Cell* left,const Cell* right)");
    Cell::init(left,right);
    Tbc_=var(*gc);
    mHbc_=var(*gc);
    Tbc_.setadata(0,gc->T0());

    vars.push_back(&Tbc_);
    vars.push_back(&mHbc_);
  }
}                // namespace tube

