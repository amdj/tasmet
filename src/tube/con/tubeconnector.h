// tubeconnector.h
//
// Author: J.A. de Jong 
//
// Description:
// This connector is able to connect two tubes together.
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef _TUBECONNECTOR_H_
#define _TUBECONNECTOR_H_
#include <array>
#include "connector.h"
#include "constants.h"


namespace tube{

  #ifndef SWIG
  class Tube;
  #endif
  #ifdef SWIG
  %catches(std::exception,...) SimpleTubeConnector::SimpleTubeConnector(us,Pos,us,Pos);
    
  #endif // SWIG
  class SimpleTubeConnector:public segment::Connector{
    std::array<us,2> segnrs;
    std::array<Pos,2> pos;
    std::array<const Tube*,2> tubes;
  public:
    SimpleTubeConnector(us seg1,Pos pos1,us seg2,Pos pos2);
    SimpleTubeConnector(const SimpleTubeConnector&);
    virtual segment::Connector* copy() const {return new SimpleTubeConnector(*this);}

    virtual vd error() const;
    #ifndef SWIG
    virtual void init(const tasystem::TaSystem& sys);

    virtual void setEqNrs(us firstdofnr);
    virtual us getNEqs() const;    
    virtual void show(us) const;
    virtual void jac(tasystem::Jacobian&) const;
    virtual void updateNf();
    #endif
  };

} // namespace tube

#endif /* _TUBECONNECTOR_H_ */
//////////////////////////////////////////////////////////////////////
