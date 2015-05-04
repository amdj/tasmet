#include "vtypes.h"
#include "constants.h"
#include "exception.h"

TRACETHIS

const char* toString(Varnr v) {
  switch(v) {
  case Varnr::none:
    return "none";
  case Varnr::rho:
    return "rho";
  case Varnr::m:
    return "m";
  case Varnr::T:
    return "T";
  case Varnr::p:
    return "p";
  case Varnr::Ts:
    return "Ts";
  case Varnr::mH:
    return "mH";
  case Varnr::U:
    return "U";
  case Varnr::Q:
    return "Q";
  case Varnr::Qs:
    return "Qs";
  default:
    throw MyError("Unknown varnr");
  }
}
