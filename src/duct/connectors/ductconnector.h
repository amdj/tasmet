// ductconnector.h
//
// Author: J.A. de Jong 
//
// Description:
// This connector is able to connect two ducts together.
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef _DUCTCONNECTOR_H_
#define _DUCTCONNECTOR_H_
#include <array>
#include "connector.h"
#include "constants.h"


namespace duct{

  #ifndef SWIG
  class Duct;
  class BcCell;
  #endif
  #ifdef SWIG
  %catches(std::exception,...) DuctConnector::DuctConnector(us,Pos,us,Pos,d K1to2=0,d K2to1=0);
  #endif // SWIG

  
  class DuctConnector:public segment::Connector{
    std::array<string,2> segids;
    std::array<Pos,2> pos;
    std::array<const BcCell*,2> bccells;
    std::array<d,2> out={{1.0f, 1.0f}};
    // Minor loss coefficients
    d K1to2, K2to1;
    // Nusselt numbers relating heat transfer at the end of a Duct
    // with solid to the temperature difference. In a later stadium,
    // these should become dependent on (at least) the Reynolds
    // number. For now they are assumed constant and equal to 1.

    // The Nusselt number is defined as
    //
    //          h * rhs
    //  Nu= --------
    //          kappa_f
    //
    // The hydraulic radius of the solid is, for the sake of
    // simplicity computed as
    // rhs= rh*Ss/Sf
    // 
    // And the heat transfer at the segment end is
    // Qs = h*Ss*(Tsbc2-Tbc1). See also documentation.
    d Nu1=1,Nu2=1;

    us firsteqnr;
  public:
    DuctConnector(const string& seg1,Pos pos1,const string& seg2,Pos pos2,d K1to2=0,d K2to1=0);
    #ifndef SWIG
    DuctConnector(const DuctConnector&)=delete;
    DuctConnector& operator=(const DuctConnector&)=delete;
    #endif // ifndef SWIG
    DuctConnector(const DuctConnector&,const tasystem::TaSystem&);    
    virtual segment::Connector* copy(const tasystem::TaSystem& s) const {return new DuctConnector(*this,s);}
    ~DuctConnector(){}
    virtual vd error() const;
    #ifndef SWIG
    virtual void setEqNrs(us firstdofnr);
    virtual us getNEqs() const;    
    virtual void show(us) const;
    virtual void jac(tasystem::Jacobian&) const;
    virtual void updateNf();
    #endif
  };

} // namespace duct

#endif /* _DUCTCONNECTOR_H_ */
//////////////////////////////////////////////////////////////////////
