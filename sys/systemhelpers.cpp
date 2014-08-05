#include "system.h"
#include "systemhelpers.h"
#include "isentropictube.h"
#include "laminarduct.h"
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
  using tube::IsentropicTube;
  using segment::connectpos;
  using segment::SegBase;
  

  

  void copyallsegsbc(TAsystem& to,const TAsystem& from){
    TRACE(14,"copyallsegsbc()");
    for(us i=0;i<from.getNSegs();i++)
      if(from[i]!=NULL)
	to.addSeg(*from[i]);
      else
	TRACE(14,"Warning! Segment is NULL");
    for(us i=0;i<from.getNBc();i++)
      if(from.getBc(i)!=NULL)
	to.addBc(*from.getBc(i));
      else
	TRACE(14,"Warning! Bc is NULL");
    
  }

  void connectbc(Seg& seg,const BcVertex& bc){
    TRACE(14,"connectbc()");
    if(bc.connectPos()==connectpos::left)
      seg.setLeftBc(vertexfrombc(copybc(bc)));
    else if(bc.connectPos()==connectpos::right)
      seg.setRightBc(vertexfrombc(copybc(bc)));
    else
      cout << "WARNING: bconnectbc(): Bc type not understood!\n";
  }

  Seg* copyseg(const Seg& orig)
  {
    if(orig.gettype().compare("IsentropicTube")==0){
      TRACE(10,"New "<<orig.gettype()<<" added to system.");
      return new IsentropicTube(static_cast<const IsentropicTube&>(orig));
    }
    else if(orig.gettype().compare("LaminarDuct")==0){
      TRACE(10,"New "<<orig.gettype()<<" added to system.");
      return new tube::LaminarDuct(static_cast<const tube::LaminarDuct&>(orig));
    }
    else if(orig.gettype().compare("LaminarDuct_e")==0){
      TRACE(10,"New "<<orig.gettype()<<" added to system.");
      return new tube::LaminarDuct_e(static_cast<const tube::LaminarDuct_e&>(orig));
    }
    else{
      TRACE(50,"Warning: segment " << orig.gettype() << " not yet implemented! Aborting...");
      abort();
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
