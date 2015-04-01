#include "bccell.h"
#include "var.h"
#include "jacrow.h"
#include "weightfactors.h"


namespace tube{
  using tasystem::JacRow;
  using tasystem::JacCol;

  BcCell::BcCell(us i,const Tube& t):
    Cell(i,t)
  {}

}                // namespace tube

