#pragma once
#ifndef _PRESSUREBC_H_
#define _PRESSUREBC_H_
#include "var.h"
#include "momentumeq.h"
#include "energyeq.h"
#include "isentropiceq.h"
#include "stateeq.h"
#include "tubebcvertex.h"

namespace tube{
  SPOILNAMESPACE
  using variable::var;

  class LeftPressureMomentumEq:public Momentum{
  public:
    virtual JacCol dpL(const TubeVertex&) const;
    virtual TubeEquation* copy() const{return new LeftPressureMomentumEq(*this);}
  };
  class LeftPressureEnergyEq:public Energy{
  public:
    virtual JacCol dpL(const TubeVertex&) const;
    virtual TubeEquation* copy() const{return new LeftPressureEnergyEq(*this);}
  };
  class LeftPressureIsentropicEq:public Isentropic{
  public:
    virtual JacCol dpL(const TubeVertex&) const;
    virtual TubeEquation* copy() const{return new LeftPressureIsentropicEq(*this);}
  };
  class LeftPressureStateEq: public State
   {
  public:
    virtual JacCol dpL(const TubeVertex&) const;
    virtual TubeEquation* copy() const{return new LeftPressureStateEq(*this);}
  };
  class LeftPressure:public TubeBcVertex
  {
    variable::var pLbc;			// Pressure boundary condition
    variable::var TLbc;			// Temperature boundary conditions
    LeftPressureMomentumEq lmomeq;
    LeftPressureEnergyEq leneq;
    LeftPressureIsentropicEq liseq;
    LeftPressureStateEq lseq;        
  public:
    // PressureBcEnergyEq peq;
    virtual void show() const;

    virtual const variable::var& pL() const final { return pLbc;}    
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
