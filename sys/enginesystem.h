#pragma once
#ifndef _ENGINESYSTEM_H_
#define _ENGINESYSTEM_H_

#include "tasystem.h"

namespace tasystem{

  class EngineSystem:public taSystem{
    virtual esdmat jac();
    virtual evd getRes();
    virtual evd setRes(evd res);
  };

} // namespace tasystem


#endif /* _ENGINESYSTEM_H_ */


