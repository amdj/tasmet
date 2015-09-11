#include "weightfactors.h"
#include "cell.h"

namespace duct{
  using std::tuple;
  using std::make_tuple;

  WeightFactors::WeightFactors(const Cell& v) {
    TRACE(15,"WeightFactors()");
    const Cell* left=v.left();
    const Cell* right=v.right();

    d vx=v.vx;
    d xl=v.xl;
    d xr=v.xr;
    // d Wfo=v.gc->getWfo();	// Wfirstorder
    if(v.left()){
      d vxim1=v.left()->vx;
      d WlLsecondorder=(vx-xl)/(vx-vxim1);
      WlL=WlLsecondorder;//-0.5*Wfo+(0.5-WlLsecondorder)*pow(Wfo,2);
      WlR=1-WlL;
    }
    else{
      WlL=1;
      WlR=0;
    }
    if(v.right()){    
      d vxip1=v.right()->vx;
      d WrRsecondorder=(xr-vx)/(vxip1-vx);
      WrR=WrRsecondorder;//+0.5*Wfo+(0.5-WrRsecondorder)*pow(Wfo,2);
      WrL=1-WrR;
    }
    else{
      WrL=0;
      WrR=1;
    }
  }
  void WeightFactors::show() const {
    TRACE(15,"WeightFactors::show()");
    cout << "WlL     :"<<WlL<<"\n";
    cout << "WlR     :"<<WlR<<"\n";
    cout << "WrL     :"<<WrL<<"\n";
    cout << "WrR     :"<<WrR<<"\n";

  }
  tuple<d,d> BcWeightFactorsW(const Cell& v){
    TRACE(15,"BcWeightfactorsW");
    assert((!v.left() && v.right()) || (v.left() && !v.right()));
    if(!v.left()){
      // Leftmost node
      d xl=0;
      d xr=v.xr;
      d xrR=v.right()->xr;

      // W0: weight factor for contribution of quantity at right
      // cell wall for something at the left cell wall

      // W1: weight factor for contribution of quantity at right
      // cell wall of neighbouring cell to what happens at the left
      // cell wall of this cell
      d W1=-xr/(xrR-xr);
      d W0=1-W1;
      // VARTRACE(40,W0);
      // VARTRACE(40,W1);
      return make_tuple(W0,W1);
    }
    else{
      d xr=v.xr;
      d xl=v.xl;
      d xlL=v.left()->xl;
      d WR2=(xr-xl)/(xlL-xl);
      d WR1=1-WR2;
      return make_tuple(WR1,WR2);
    }
  }
  tuple<d,d> BcWeightFactorsV(const Cell& v){
    TRACE(15,"BcWeightFactorsV()");
    assert((!v.left() && v.right()) || (v.left() && !v.right()));
    if(!v.left()){
      // Leftmost node
      d xl=0;
      d xr=v.vx;
      d xrR=v.right()->vx;

      d W1=-xr/(xrR-xr);
      d W0=1-W1;
      return make_tuple(W0,W1);
    }
    else{
      d xr=v.xr;
      d xl=v.vx;
      d xlL=v.left()->vx;
      d WR2=(xr-xl)/(xlL-xl);
      d WR1=1-WR2;
      return make_tuple(WR1,WR2);
    }
  }

} // namespace duct



