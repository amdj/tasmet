#include "segconnection.h"
#include "system.h"

namespace tasystem{

  void coupleSegs(const SegConnection& sc,TAsystem& sys){
    TRACE(14,"coupleSegs()");
    Seg& seg1=*sys.getSeg(sc.firstseg);
    Seg& seg2=*sys.getSeg(sc.secondseg);
    SegCoupling coupling=sc.coupling;
    us seg1size=seg1.vvertex.size();
    us seg2size=seg2.vvertex.size();
    assert(seg1size>0);
    assert(seg2size>0);
    if (coupling==tailhead){
      // Seg1 is coupled with its tail to Seg2's head
      TRACE(3,"Coupling seg1 with its tail to the head of seg2");
      seg1.setRight(seg2);
      seg2.setLeft(seg1);
      seg1.vvertex[seg1size-1]->right=seg2.vvertex[0].get();
      seg2.vvertex[0]->left=seg1.vvertex[seg1size-1].get();

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
