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
  class BcCell;
  #endif
  #ifdef SWIG
  %catches(std::exception,...) SimpleTubeConnector::SimpleTubeConnector(us,Pos,us,Pos,d K1to2=0,d K2to1=0);
  #endif // SWIG

  
  class SimpleTubeConnector:public segment::Connector{
    std::array<string,2> segids;
    std::array<Pos,2> pos;
    std::array<const BcCell*,2> bccells;
    std::array<d,2> out={{1.0f, 1.0f}};
    // Minor loss coefficients
    d K1to2, K2to1;
    us firsteqnr;
  public:
    SimpleTubeConnector(const string& seg1,Pos pos1,const string& seg2,Pos pos2,d K1to2=0,d K2to1=0);
    SimpleTubeConnector(const SimpleTubeConnector&)=delete;
    SimpleTubeConnector(const SimpleTubeConnector&,const tasystem::TaSystem&);    
    virtual segment::Connector* copy(const tasystem::TaSystem& s) const {return new SimpleTubeConnector(*this,s);}
    ~SimpleTubeConnector(){}
    virtual vd error() const;
    #ifndef SWIG
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
