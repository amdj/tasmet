// laminarduct.h
//
// Author: J.A. de Jong 
//
// Description:
// 
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef LAMINARDUCT_H_
#define LAMINARDUCT_H__
#include "tube.h"
#include "vtypes.h"
#include "laminardrag.h"
#include "hopkinsheat.h"

namespace solids{
  class Solid;
}

namespace tube{
  #ifndef SWIG
  SPOILNAMESPACE
  #endif
  #ifdef SWIG
  %catches(std::exception,...) LaminarDuct::LaminarDuct(const Geom&,const string& solid);
  %catches(std::exception,...) LaminarDuct::LaminarDuct(const Geom&);
  %catches(std::exception,...) LaminarDuct::setSolid(const string& solid);
  #endif

  class LaminarDuct:public Tube
  {
    drag::LaminarDragResistance laminardrag;
    HopkinsHeatSource hopkinsheat;
    // If isolated, no time-averaged heat transfer is allowed between
    // tube and
    bool insulated=false;
    // If present a solid
    const solids::Solid* solid=nullptr;

    // Temporary for testing
    d Tl=-1,Tr=-1;

  protected:
    LaminarDuct(const LaminarDuct&,const tasystem::TaSystem&);
  public:
    LaminarDuct(const Geom& geom,d Tl,d Tr);
    LaminarDuct(const Geom& geom);
    LaminarDuct(const Geom& geom,const string& solid);
    LaminarDuct(const LaminarDuct&)=delete;
    void init();
    void setInsulated(bool i){insulated=i;}
    bool isInsulated() const {return insulated;}
    void setSolid(const string& solid,d ksfrac=1);
    bool hasSolid() const {return solid?true:false;}
    const solids::Solid& getSolid() const;
    void setVarsEqs(Cell&) const;
    LaminarDuct& operator=(const LaminarDuct&)=delete;
    void show(us) const;
    segment::Seg* copy(const tasystem::TaSystem& s) const {return new LaminarDuct(*this,s);}
    #ifndef SWIG
    virtual ~LaminarDuct();
    const HeatSource& heatSource() const { return hopkinsheat;}
    virtual const drag::DragResistance& dragResistance() const {return laminardrag;}
    #endif
  };

  
} /* namespace tube */

#endif	// LAMINARDUCT_H_






