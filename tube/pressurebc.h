#pragma once
#ifndef _PRESSUREBC_H_
#define _PRESSUREBC_H_
#include "bcvertex.h"

namespace tube{

  class LeftPressure:public TubeBcVertex
  {
  public:
    LeftPressure(const Tube& t);
    LeftPressure(const Tube&t,variable::var& pres);
    LeftPressure(const Tube&t,variable::var& pres,variable::var& temp);

    ~LeftPressure(){}
    void updateW();

    virtual vd msource() const;	// Prescribed left pressure
    virtual vd esource() const;	// Same
    variable::var pL;			// Pressure boundary condition
    variable::var TL;			// Temperature boundary conditions
  };

} // namespace tube

#endif /* _PRESSUREBC_H_ */
