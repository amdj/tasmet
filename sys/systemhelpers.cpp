#include "systemhelpers.h"
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
      seg.setLeftbc(vertexfrombc(copybc(bc)));
    else if(bc.connectPos()==connectpos::right)
      seg.setRightbc(vertexfrombc(copybc(bc)));
    else
      cout << "WARNING: bconnectbc(): Bc type not understood!\n";
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
    TRACE(14,"copbybc(BcVertex&)");
    if(orig.gettype().compare("RightImpedance")==0){
      TRACE(14,"Copying a RightImpedance bc...");
      return static_cast<BcVertex*>(new RightImpedance(static_cast<const RightImpedance&>(orig)));
    } else if(orig.gettype().compare("LeftPressure")==0){
      TRACE(14,"Copying a LeftPressure bc...");
      return static_cast<BcVertex*>(new LeftPressure(static_cast<const LeftPressure&>(orig)));
    }
    else{
      TRACE(50,"Boundary condition type not understood");
      return NULL;
      exit(1);
    }
  }
  Vertex* vertexfrombc(BcVertex* orig){
    if(orig->gettype().compare("RightImpedance")==0)
      return static_cast<Vertex*>(static_cast<RightImpedance*>(orig));
    else if(orig->gettype().compare("LeftPressure")==0)
      return static_cast<Vertex*>(static_cast<LeftPressure*>(orig));
    else{
      TRACE(50,"Boundary condition type not understood");
      return NULL;
      exit(1);
    }
  }
  
  
}		// namespace tasystem
