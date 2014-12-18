#include "seg.h"
#include "tasystem.h"
#include <cassert>
namespace segment{

  Seg::Seg(const Seg& other):SegConBase(other){
    TRACE(13,"Seg::Seg()");
  }
  Seg::Seg():SegConBase(){}

}		 // Namespace segment
