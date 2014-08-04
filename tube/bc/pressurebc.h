#pragma once
#ifndef _PRESSUREBC_H_
#define _PRESSUREBC_H_

#include "energyeq.h"
#include "tubebcvertex.h"

namespace tube{

  using variable::var;
  using segment::connectpos;

  class PressureBcEnergyEq:public Energy
  {
  public:
    virtual dmat dUi(const TubeVertex& v) const;
    virtual dmat dUip1(const TubeVertex& v) const;
  };
  
  class LeftPressure:public TubeBcVertex
  {
  public:
    variable::var pL;			// Pressure boundary condition
    variable::var TL;			// Temperature boundary conditions

    PressureBcEnergyEq peq;    
    LeftPressure(us,const var&);
    LeftPressure(us,const var&,const var& temp);
    LeftPressure(const LeftPressure& other);
    LeftPressure& operator=(const LeftPressure&);
    ~LeftPressure(){}
    virtual void initTubeVertex(us i,const Tube&);
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
