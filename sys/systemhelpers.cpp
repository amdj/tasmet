#include "systemhelpers.h"
#include "bcvertex.h"
#include "tube.h"
#include "impedancebc.h"
#include "pressurebc.h"

namespace tasystem{
  using tube::Tube;
  using tube::RightImpedance;
  using tube::LeftPressure;
  using segment::connectpos;

  void connectbc(Seg& seg,const BcVertex& bc){
    TRACE(14,"connectbc()");
    if(bc.connectPos()==connectpos::left)
      seg.setLeftbc((Vertex&) bc);
    else if(bc.connectPos()==connectpos::right)
      seg.setRightbc((Vertex&) bc);
    else
      cout << "WARNING: connectbc(): Bc type not understood!\n";
  }

  Seg* copyseg(const Seg& orig)
  {
    if(orig.gettype().compare("Tube")==0){
      TRACE(10,"New tube added to system.");
      return new Tube((const Tube&) orig);
    }
    else{
      TRACE(50,"Warning: segment " << orig.gettype() << " not yet implemented!");
      return NULL;
    }
    
  } // copyseg
  BcVertex* copybc(const BcVertex& orig){
    if(orig.gettype().compare("RightImpedance")==0){
      return new RightImpedance((RightImpedance&) orig);
    }
    else if(orig.gettype().compare("LeftPressure")==0){
      return new LeftPressure((LeftPressure&) orig);
    }

    else{
      TRACE(50,"Boundary condition type not understood");
      return NULL;
    }

  }
  
  
}		// namespace tasystem
