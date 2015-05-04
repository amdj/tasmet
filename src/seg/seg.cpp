#include "seg.h"
#include "tasystem.h"
#include <cassert>
#include "exception.h"
namespace segment{
  using tasystem::TaSystem;

  Seg::Seg():SegConBase(){}
  void Seg::setPhaseContraint(tasystem::PhaseConstraint) {
    TRACE(15,"tasystem::PhaseConstraint v)");
    throw MyError("This segment is unable to constrain the phase on.");
  }
}		 // Namespace segment
