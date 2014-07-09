#pragma once
#ifndef _PRESSUREBC_H_
#define _PRESSUREBC_H_


#include "tubebcvertex.h"

namespace tube{

  using tasystem::Globalconfptr;  
  using variable::var;

  using segment::connectpos;
  
  class LeftPressure:public TubeBcVertex
  {
  public:
    variable::var pL;			// Pressure boundary condition
    variable::var TL;			// Temperature boundary conditions

    LeftPressure(us,var&);
    LeftPressure(us,var&,var& temp);

    ~LeftPressure(){}
    void updateW(const Geom& geom);
    virtual string gettype() const {return string("LeftPressure");}
    virtual enum connectpos connectPos() const {return connectpos::left;}
    virtual vd msource() const;	// Prescribed left pressure
    virtual vd esource() const;	// Same
  };

} // namespace tube

#endif /* _PRESSUREBC_H_ */
