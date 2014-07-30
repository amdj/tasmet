#include "system.h"
#include "systemhelpers.h"
#include "tube.h"
#include "impedancebc.h"
#include "twimpedance.h"
#include "pressurebc.h"
#include "isotwall.h"
namespace tasystem{
  using tube::Tube;
  using tube::RightImpedance;
  using tube::RightIsoTWall;  
  using tube::TwImpedance;
  using tube::LeftPressure;
  using segment::connectpos;
  using segment::SegBase;
  

  

  void copyallsegsbc(TAsystem& to,const TAsystem& from){
    TRACE(14,"copyallsegsbc()");
    for(us i=0;i<from.getNsegs();i++)
      if(from[i]!=NULL)
	to.addseg(*from[i]);
      else
	TRACE(14,"Warning! Segment is NULL");
    for(us i=0;i<from.getNbc();i++)
      if(from.getBc(i)!=NULL)
	to.addbc(*from.getBc(i));
      else
	TRACE(14,"Warning! Bc is NULL");
    
  }

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
      return new Tube(static_cast<const Tube&>(orig));
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
    }
    else if(orig.gettype().compare("TwImpedance")==0){
      TRACE(14,"Copying a TwImpedance bc...");
      return static_cast<BcVertex*>(new TwImpedance(static_cast<const TwImpedance&>(orig)));
    }
    else if(orig.gettype().compare("RightIsoTWall")==0){
      TRACE(14,"Copying a RightIsoTWall bc...");
      return static_cast<BcVertex*>(new RightIsoTWall(static_cast<const RightIsoTWall&>(orig)));
    }
    else if(orig.gettype().compare("LeftPressure")==0){
      TRACE(14,"Copying a LeftPressure bc...");
      return static_cast<BcVertex*>(new LeftPressure(static_cast<const LeftPressure&>(orig)));
    }
    else{
      WARN("Boundary condition type not understood");
      abort();
      return NULL;		// For the compiler. But we quit executing.
    }
  }
  Vertex* vertexfrombc(BcVertex* orig){
    if(orig->gettype().compare("RightImpedance")==0)
      return static_cast<Vertex*>(static_cast<RightImpedance*>(orig));
    else if(orig->gettype().compare("TwImpedance")==0)
      return static_cast<Vertex*>(static_cast<TwImpedance*>(orig));
    else if(orig->gettype().compare("RightIsoTWall")==0)
      return static_cast<Vertex*>(static_cast<RightIsoTWall*>(orig));
    else if(orig->gettype().compare("LeftPressure")==0)
      return static_cast<Vertex*>(static_cast<LeftPressure*>(orig));
    else{
      WARN("Boundary condition type not understood");
      abort();
    }
  }
  
  
}		// namespace tasystem
