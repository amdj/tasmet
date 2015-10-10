// ductpistonconnector.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef DUCTPISTONCONNECTOR_H
#define DUCTPISTONCONNECTOR_H

#include <array>
#include "connector.h"
#include "constants.h"


namespace duct {
  
  #ifndef SWIG
  class Duct;
  class BcCell;
  #endif
  
  
} // namespace duct


namespace mech{
  #ifndef SWIG
  class Piston;
  #endif // SWIG  

  #ifdef SWIG
  %catches(std::exception,...) DuctPistonConnector::DuctPistonConnector(us ductnr,Pos ductpos,us pistonnr,Pos pistonpos,d KDuctPiston=0,d KPistonDuct=0);
  #endif
  
  class DuctPistonConnector:public segment::Connector{
    string pistonid,ductid;
    Pos pistonPos,ductPos;
    us firsteqnr;
    // d dx=0,Sfgem;
    const duct::Duct* duct=nullptr;
    const Piston* piston=nullptr;

    d KDuctPiston;              // Minor loss coefficient for flow
                                // from the duct to the piston volume
    d KPistonDuct;              // Minor loss coefficient for flow
                                // from the piston volume to the duct
  
    // See documentation for more information.
    
  public:
    DuctPistonConnector(const string& ductid,Pos ductpos,const string& pistonid,Pos pistonpos,d KDuctPiston=0,d KPistonDuct=0);
    #ifndef SWIG
    DuctPistonConnector(const DuctPistonConnector&)=delete;
    DuctPistonConnector& operator=(const DuctPistonConnector&)=delete;
    #endif // ifndef SWIG
    DuctPistonConnector(const DuctPistonConnector&,const tasystem::TaSystem&);    
    virtual segment::Connector* copy(const tasystem::TaSystem& s) const {return new DuctPistonConnector(*this,s);}
    ~DuctPistonConnector(){}
    virtual vd error() const;
    #ifndef SWIG
    virtual void setEqNrs(us firstdofnr);
    virtual us getNEqs() const;    
    virtual void show(us) const;
    virtual void jac(tasystem::Jacobian&) const;
    virtual void updateNf();
  private:
    vd kappaSft() const;
    #endif
  };

} // namespace mech

#endif /* _DUCTCONNECTOR_H_ */
//////////////////////////////////////////////////////////////////////
