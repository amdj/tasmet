#include "weightfactors.h"
#include "cell.h"

namespace tube{

  std::tuple<d,d,d,d> WeightFactors(const Cell& v) {
    TRACE(15,"WeightFactors()");
    const Cell* left=v.left();
    const Cell* right=v.right();

    d vx=v.vx;
    d xL=v.xL;
    d xR=v.xR;

    d WRr=0,WRl=0,WLr=0,WLl=0;
    
    if(v.left()){
      d vxim1=v.left()->vx;
      WLl=(vx-xL)/(vx-vxim1);
      WLr=1-WLl;
    }
    else{
      WLl=1;
      WLr=0;
    }
    if(v.right()){    
      d vxip1=v.right()->vx;
      WRr=(xR-vx)/(vxip1-vx);
      WRl=1-WRr;
    }
    else{
      WRl=0;
      WRr=1;
    }
    return std::make_tuple(WLl,WLr,WRl,WRr);
  }
} // namespace tube



