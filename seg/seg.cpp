#include "seg.h"

namespace segment{
  using tasystem::Globalconf;   

  Seg::Seg(){
    TRACE(13,"Seg::Seg()");
  }
  Seg::Seg(const Seg& other){}
  void Seg::init(const Globalconf& gc1){
    TRACE(13,"Seg::init()");
    this->gc=&gc1;
    init_=true;
  }  
}		 // Namespace segment
