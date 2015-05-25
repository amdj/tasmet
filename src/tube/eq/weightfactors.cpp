#include "weightfactors.h"
#include "cell.h"

namespace tube{
  using std::tuple;
  using std::make_tuple;

  WeightFactors::WeightFactors(const Cell& v) {
    TRACE(15,"WeightFactors()");
    const Cell* left=v.left();
    const Cell* right=v.right();

    d vx=v.vx;
    d xL=v.xL;
    d xR=v.xR;

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
  }
  tuple<d,d> BcWeightFactorsW(const Cell& v){
    TRACE(15,"BcWeightfactorsW");
    assert((!v.left() && v.right()) || (v.left() && !v.right()));
    if(!v.left()){
      // Leftmost node
      d xL=0;
      d xR=v.xR;
      d xRR=v.right()->xR;

      // W0: weight factor for contribution of quantity at right
      // cell wall for something at the left cell wall

      // W1: weight factor for contribution of quantity at right
      // cell wall of neighbouring cell to what happens at the left
      // cell wall of this cell
      d W1=-xR/(xRR-xR);
      d W0=1-W1;
      // VARTRACE(40,W0);
      // VARTRACE(40,W1);
      return make_tuple(W0,W1);
    }
    else{
      d xR=v.xR;
      d xL=v.xL;
      d xLL=v.left()->xL;
      d WR2=(xR-xL)/(xLL-xL);
      d WR1=1-WR2;
      return make_tuple(WR1,WR2);
    }
  }
  tuple<d,d> BcWeightFactorsV(const Cell& v){
    TRACE(15,"BcWeightFactorsV()");
    assert((!v.left() && v.right()) || (v.left() && !v.right()));
    if(!v.left()){
      // Leftmost node
      d xL=0;
      d xR=v.vx;
      d xRR=v.right()->vx;

      d W1=-xR/(xRR-xR);
      d W0=1-W1;
      return make_tuple(W0,W1);
    }
    else{
      d xR=v.xR;
      d xL=v.vx;
      d xLL=v.left()->vx;
      d WR2=(xR-xL)/(xLL-xL);
      d WR1=1-WR2;
      return make_tuple(WR1,WR2);
    }
  }

} // namespace tube



