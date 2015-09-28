// laminarduct.h
//
// Author: J.A. de Jong 
//
// Description:
// LaminarDuct derives 'virtual public' from Duct, because a stack
// contains both a DuctWithSolid and a LaminarDuct. See also documentation
// on DuctWithSolid.
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef LAMINARDUCT_H_
#define LAMINARDUCT_H__
#include "duct.h"
#include "vtypes.h"

namespace duct{
  #ifndef SWIG
  SPOILNAMESPACE
  #endif
  #ifdef SWIG
  %catches(std::exception,...) LaminarDuct::LaminarDuct(const Geom&);
  %catches(std::exception,...) LaminarDuct::LaminarDuct(const Geom& geom,d TwL,d Twr=-1);
  %catches(std::exception,...) LaminarDuct::setInsulated(bool i);
  #endif

  #ifndef SWIG
  namespace drag {
    class DragResistance;
  } // namespace drag
  #endif

  class LaminarDuct:virtual public Duct
  {
    drag::DragResistance* laminardrag=nullptr;
    HeatSource* hopkinsheat=nullptr;
    // If isolated, no time-averaged heat transfer is allowed between
    // duct and
    bool insulated=false;

    // Initialize temperatures with this value, when nonzero
    vd Tinit;

  protected:
    LaminarDuct(const LaminarDuct&,const tasystem::TaSystem&);
  public:
    LaminarDuct(const Geom& geom,d TwL,d Twr=-1);
    LaminarDuct(const Geom& geom);
    LaminarDuct(const LaminarDuct&)=delete;
    LaminarDuct& operator=(const LaminarDuct&)=delete;    
    virtual ~LaminarDuct();

    // Set and get the insulated flag.
    void setInsulated(bool i);
    bool isInsulated() const {return insulated;}
    segment::Seg* copy(const tasystem::TaSystem& s) const {return new LaminarDuct(*this,s);}
    
    #ifndef SWIG

    void init();
    void show(us) const;

    void setVarsEqs(Cell& c) const;
    const HeatSource& heatSource() const {return *hopkinsheat;}
    const drag::DragResistance& dragResistance() const {return *laminardrag;}
    #endif

  };

  
} /* namespace duct */

#endif	// LAMINARDUCT_H_






