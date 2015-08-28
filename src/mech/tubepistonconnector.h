// tubepistonconnector.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef TUBEPISTONCONNECTOR_H
#define TUBEPISTONCONNECTOR_H

#include <array>
#include "connector.h"
#include "constants.h"


namespace tube {
  
  #ifndef SWIG
  class Tube;
  class BcCell;
  #endif
  
  
} // namespace tube


namespace mech{
  #ifndef SWIG
  class Piston;
  #endif // SWIG  

  #ifdef SWIG
  %catches(std::exception,...) TubePistonConnector::TubePistonConnector(us tubenr,Pos tubepos,us pistonnr,Pos pistonpos,d KTubePiston=0,d KPistonTube=0);
  #endif
  
  class TubePistonConnector:public segment::Connector{
    string pistonid,tubeid;
    Pos pistonPos,tubePos;
    us firsteqnr;
    // d dx=0,Sfgem;
    const tube::Tube* tube=nullptr;
    const Piston* piston=nullptr;

    d KTubePiston;              // Minor loss coefficient for flow
                                // from the tube to the piston volume
    d KPistonTube;              // Minor loss coefficient for flow
                                // from the piston volume to the tube
  
    // See documentation for more information.
    
  public:
    TubePistonConnector(const string& tubeid,Pos tubepos,const string& pistonid,Pos pistonpos,d KTubePiston=0,d KPistonTube=0);
    TubePistonConnector(const TubePistonConnector&)=delete;
    TubePistonConnector(const TubePistonConnector&,const tasystem::TaSystem&);    
    virtual segment::Connector* copy(const tasystem::TaSystem& s) const {return new TubePistonConnector(*this,s);}
    ~TubePistonConnector(){}
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

#endif /* _TUBECONNECTOR_H_ */
//////////////////////////////////////////////////////////////////////
