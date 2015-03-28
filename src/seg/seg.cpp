#include "seg.h"
#include "tasystem.h"
#include <cassert>
namespace segment{
  using tasystem::TaSystem;

  Seg::Seg(const Seg& other):SegConBase(other){
    TRACE(13,"Seg::Seg()");
  }
  void Seg::init(const TaSystem& sys){
    TRACE(15,"Seg::init()");
    SegConBase::init(sys);
  }
  Seg::Seg():SegConBase(){}

}		 // Namespace segment
