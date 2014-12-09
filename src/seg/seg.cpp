#include "seg.h"
#include "tasystem.h"
#include <cassert>
namespace segment{
  using tasystem::Globalconf;   
  using tasystem::TaSystem;   

  Seg::Seg(){
    TRACE(13,"Seg::Seg()");
  }
  Seg::Seg(const Seg& other){}
  void Seg::init(const TaSystem& sys){
    TRACE(13,"Seg::init()");
    this->gc=&sys.gc;
    init_=true;
  }  
  const Globalconf& Seg::getGc() const{
    assert(gc);
    return *gc;
  }
}		 // Namespace segment
