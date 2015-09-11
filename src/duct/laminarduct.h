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
#include "duct.h"
#include "vtypes.h"
#include "laminardrag.h"
#include "hopkinsheat.h"

namespace solids{
  class Solid;
}

namespace duct{
  #ifndef SWIG
  SPOILNAMESPACE
  #endif
  #ifdef SWIG
  %catches(std::exception,...) LaminarDuct::LaminarDuct(const Geom&,const string& solid);
  %catches(std::exception,...) LaminarDuct::LaminarDuct(const Geom&);
  %catches(std::exception,...) LaminarDuct::setQsin(d Qsin);
  %catches(std::exception,...) LaminarDuct::setSolid(const string& solid,d ksfrac=1);
  #endif

  class LaminarDuct:public Duct
  {
    drag::LaminarDragResistance laminardrag;
    HopkinsHeatSource hopkinsheat;
    // If isolated, no time-averaged heat transfer is allowed between
    // duct and
    bool insulated=false;
    // If present a solid
    const solids::Solid* solid=nullptr;
      // Temporary for testing
    d Tl=-1,Tr=-1;
    d Qsin=0;			// Heat input into the solid
  protected:
    LaminarDuct(const LaminarDuct&,const tasystem::TaSystem&);
  public:
    LaminarDuct(const Geom& geom,d Tl,d Tr);
    LaminarDuct(const Geom& geom);
    LaminarDuct(const Geom& geom,const string& solid);
    LaminarDuct(const LaminarDuct&)=delete;
    void init();
    void setQsin(d Qsin);	// Apply solid heat input. For this, a
				// solid need to be defined with setSolid().
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

  
} /* namespace duct */

#endif	// LAMINARDUCT_H_






