#include "seg.h"
#include "tasystem.h"
#include <cassert>
namespace segment{

  Seg::Seg(const Seg& other):Connector(other){
    TRACE(13,"Seg::Seg()");
  }
  Seg::Seg():Connector(){}

}		 // Namespace segment
