// isentropictube.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef ISENTROPICTUBE_H
#define ISENTROPICTUBE_H

#include "vtypes.h"
#include "duct.h"


namespace duct{

  class IsentropicTube:public Duct
  {
    drag::DragResistance nodrag;
    HeatSource noheat;
    IsentropicTube& operator=(const IsentropicTube&) =delete;    
    IsentropicTube(const IsentropicTube&)=delete;
    IsentropicTube(const IsentropicTube&,const tasystem::TaSystem&);
  public:
    IsentropicTube(const Geom& geom);
    virtual segment::Seg* copy(const tasystem::TaSystem&) const;
    virtual ~IsentropicTube();

    // Push the right variables and equations
    virtual void setVarsEqs(Cell&) const;
    #ifndef SWIG

    vd dragCoefVec(us i) const;
    virtual const drag::DragResistance& dragResistance() const {return nodrag;}
    virtual const HeatSource& heatSource() const {return noheat;}
    void cleanup();
    #endif  // ifndef SWIG
  };
  
} /* namespace duct */


#endif // ISENTROPICTUBE_H
//////////////////////////////////////////////////////////////////////
