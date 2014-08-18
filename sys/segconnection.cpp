#include "segconnection.h"
#include "system.h"

namespace tasystem{

  void coupleSegs(const SegConnection& sc,taSystem& sys){
    TRACE(14,"coupleSegs()");
    SegBase& seg1=*sys.getSeg(sc.firstseg);
    SegBase& seg2=*sys.getSeg(sc.secondseg);
    SegCoupling coupling=sc.coupling;
    if (coupling==tailhead){
      // Seg1 is coupled with its tail to Seg2's head
      TRACE(3,"Coupling seg1 with its tail to the head of seg2");
      seg1.setRight(seg2);
      seg2.setLeft(seg1);

    }
    else if(coupling==headtail){
      // Seg2 is coupled with its tail to Seg1's head
      TRACE(3,"Coupling seg1 with its head to the tail of seg2");
      seg1.setLeft(seg2);
      seg2.setRight(seg1);
    }
    else if(coupling==tailtail){
      seg1.setRight(seg2);
      seg2.setRight(seg1);
    }
    else {
      // Coupling is headhead
      seg1.setLeft(seg2);
      seg2.setRight(seg1);
      WARN("Segment coupling specified not (yet) implemented!");
    }
  } // coupleSegs()

} // namespace tasystem
