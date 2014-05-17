#pragma once
#include <logger.h>
#include <material.h>
#include "var.h"

namespace gases{
  SPOILNAMESPACE
  using variable::var;
  class Gasvar:public Gas{
  public:
    var rho(const var& T,const var& p);
    var p(const var& T,const var& rho);
    var cp(const var& T);
    var pr(const var& T) ;
    var h(const var& T) ;
    var cv(const var& T);
    var e(const var& T) ;
    var beta(const var& T);
    var gamma(const var& T) ;
    var cm(const var& T);
    var mu(const var& T) ;
    var kappa(const var& T) ;

  };				// class Gasvar

} // Namespace gases
