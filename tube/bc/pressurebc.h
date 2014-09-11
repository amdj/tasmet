#pragma once
#ifndef _PRESSUREBC_H_
#define _PRESSUREBC_H_
#include "var.h"

#include "energyeq.h"
#include "tubebcvertex.h"

namespace tube{
  SPOILNAMESPACE
  using variable::var;

  // class PressureBcEnergyEq:public Energy
  // {
  // public:
  //   virtual dmat dUi(const TubeVertex& v) const;
  //   virtual dmat dUip1(const TubeVertex& v) const;
  // };

  class LeftPressureEq:public TubeEquation{
    var& pLbc;
    JacCol dpL(const TubeVertex&) const;
    
  public:
    LeftPressureEq(variable::var& pLbc): pLbc(pLbc){}
    virtual vd error(const TubeVertex&) const;
    virtual JacRow jac(const TubeVertex& v) const;
    virtual TubeEquation* copy() const{return new LeftPressureEq(*this);}
    
  };
  
  class LeftPressure:public TubeBcVertex
  {
    variable::var pLbc;			// Pressure boundary condition
    variable::var TLbc;			// Temperature boundary conditions
    LeftPressureEq leq;
  public:
    // PressureBcEnergyEq peq;
    virtual void show() const;

    
    LeftPressure(const var&);
    LeftPressure(const var&,const var& temp);
    LeftPressure(const LeftPressure& other);
    LeftPressure& operator=(const LeftPressure&);
    ~LeftPressure(){}
    virtual void initTubeVertex(us i,const Tube&);
    virtual TubeBcVertex* copy() const { return new LeftPressure(*this);}
  private:
    void updateW(const SegBase&);
  public:
    virtual string getType() const {return string("LeftPressure");}
    virtual enum connectpos connectPos() const {return connectpos::left;}
    virtual vd msource() const;	// Prescribed left pressure
    virtual vd esource() const;	// Same
  };

} // namespace tube

#endif /* _PRESSUREBC_H_ */
