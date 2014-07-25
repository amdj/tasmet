#pragma once
#ifndef _PRESSUREBC_H_
#define _PRESSUREBC_H_


#include "tubebcvertex.h"

namespace tube{

  using variable::var;
  using segment::connectpos;
  
  class LeftPressure:public TubeBcVertex
  {
  public:
    variable::var pL;			// Pressure boundary condition
    variable::var TL;			// Temperature boundary conditions

    LeftPressure(us,const var&);
    LeftPressure(us,const var&,const var& temp);
    LeftPressure(const LeftPressure& other);
    LeftPressure& operator=(const LeftPressure&);
    ~LeftPressure(){}
    virtual void Init(us i,const SegBase&);
  private:
    void updateW(const SegBase&);
  public:
    virtual string gettype() const {return string("LeftPressure");}
    virtual enum connectpos connectPos() const {return connectpos::left;}
    virtual vd msource() const;	// Prescribed left pressure
    virtual vd esource() const;	// Same
  };

} // namespace tube

#endif /* _PRESSUREBC_H_ */
