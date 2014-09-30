#pragma once
#ifndef _ENGINEPRESSURE_H_
#define _ENGINEPRESSURE_H_


#include "adiabaticwall.h"
#include "continuityeq.h"
#include "momentumeq.h"
#include "energyeq.h"


namespace tube{

  class LeftEnginePressure;
  template <typename T>
  class LeftEnginePressureEq:public TubeEquation {
    LeftEnginePressure* lep;
    T& N;                       // Normal equation of this type
    T& A;                       // Amplitute-prescribed version
  public:

    LeftEnginePressureEq(T& N,T& A):N(N),A(A){}
    virtual ~LeftEnginePressureEq(){}
    virtual vd error(const TubeVertex& v) const{
      TRACE(15,"LeftEnginePressureEq::error()");
      vd error=N.error(v);
      if(v.gc->Nf>0){
        TRACE(15,"Filling amplitude errror..");
        vd Aerror=A.error(v);
        error(1)=Aerror(1);
        error(2)=Aerror(2);
      }
      return error;
    }
    virtual JacRow jac(const TubeVertex& v) const{
      TRACE(15,"LeftEnginePressureEq::jac()");
      
      JacRow cjac=N.jac(v);
      cjac.setRowDof(dofnr);
      if(v.gc->Nf>0){
        us njaccols=cjac.jaccols.size();
        JacRow cAjac=A.jac(v);
        TRACE(15,"Filling amplitude errror..");
        // Now we mix up the second and third row of each Jacobian column
        assert(cjac.jaccols.size()==cAjac.jaccols.size());

        for(us i=0;i<njaccols;i++) {
          JacCol& cjacc=cjac.jaccols.at(i);
          JacCol& cAjacc=cAjac.jaccols.at(i);      
        
          cjacc.data().row(1)=cAjacc.data().row(1);
          cjacc.data().row(2)=cAjacc.data().row(2);        
        } // for
      }
      return cjac;  
    }
  };

  class LeftEnginePressureState:public TubeEquation{
    const LeftEnginePressure* lep;
  public:
    LeftEnginePressureState(const LeftEnginePressure* l):
      lep(l){}
    virtual vd error(const TubeVertex&) const;
    virtual JacRow jac(const TubeVertex& v) const;
  private:
  private:
    JacCol dpL(const TubeVertex&) const;
    JacCol drhoi(const TubeVertex&) const;
    JacCol dTi(const TubeVertex&) const;
    JacCol drhoim1(const TubeVertex&) const;
    JacCol dTim1(const TubeVertex&) const;
    JacCol drhoip1(const TubeVertex&) const;
    JacCol dTip1(const TubeVertex&) const;
  };
  
  
  
  class LeftEnginePressure:public LeftAdiabaticWall{
    d amplitude__;
    Continuity cA;    
    LeftEnginePressureEq<Continuity> lepc;
    Momentum mA;    
    LeftEnginePressureEq<Momentum> lepm;
    Energy eA;    
    LeftEnginePressureEq<Energy> lepe;
    LeftEnginePressureState sA;
    LeftEnginePressureEq<TubeEquation> leps;    


  public:
    LeftEnginePressure(d a):
      amplitude__(a),
      lepc(c,cA),
      lepm(m,mA),
      lepe(e,eA),
      sA(this),
      leps(sL,sA){}
    LeftEnginePressure(const LeftEnginePressure& o): LeftEnginePressure(o.amplitude__){}
    d amplitude() const {return amplitude__;}
    void setAmplitude(d a) {amplitude__=a;}
    virtual TubeBcVertex* copy() const {return new LeftEnginePressure(*this);}
    virtual enum connectpos connectPos() const {return connectpos::left;}
    virtual void initTubeVertex(us i,const Tube& thisseg);

  };


} // namespace tube

#endif /* _ENGINEPRESSURE_H_ */
