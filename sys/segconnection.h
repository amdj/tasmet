#pragma once
#ifndef _SEGCONNECTION_H_
#define _SEGCONNECTION_H_
#include <vtypes.h>




namespace tasystem{
SPOILNAMESPACE
  enum SegCoupling{
    headhead,tailtail,headtail,tailhead
  };
  class TAsystem;
  class SegConnection
  {
  public:
    SegConnection(us seg1,us seg2,SegCoupling s): firstseg(seg1),secondseg(seg2),coupling(s) {}
    us firstseg=0;
    us secondseg=0;
    SegCoupling coupling;
    virtual ~SegConnection() {}
  };  
  
  void coupleSegs(const SegConnection&,TAsystem& sys); // Couple two segments  


} // namespace tasystem



#endif /* _SEGCONNECTION_H_ */
