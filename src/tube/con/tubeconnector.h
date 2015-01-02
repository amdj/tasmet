#pragma once
#ifndef _TUBECONNECTOR_H_
#define _TUBECONNECTOR_H_

#include "connector.h"


namespace tube{
  // enum connectpos{ left,right};	// Where to connect the boundary condition.

  class TubeConnector:public segment::Connector{
  protected:
    vector<us> segnrs;
    vector<position> connectpos;
    vector<const Tube*> tubes;
  public:
    TubeConnector()
  };
  
  class SimpleTubeConnector:public TubeConnector{
  public:
    SimpleTubeConnector(us seg1,us seg2,position pos1,position pos2);
  };

} // namespace tube

#endif /* _TUBECONNECTOR_H_ */
