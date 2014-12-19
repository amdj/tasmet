#include "weightfactors.h"
#include "tubevertex.h"
#include "tube.h"

namespace tube{


  WeightFactors::WeightFactors(const TubeVertex& v):
    LocalGeom(v.getTube().geom(),v.geti())
  {
    TRACE(15,"WeightFactors::WeightFactors()");
    const TubeVertex* left=v.left();
    const TubeVertex* right=v.right();

    if(left) {   
      TRACE(15,"left present");
      const LocalGeom llg=LocalGeom(*left);
      vSfL=llg.vSf;
      vxm1=llg.vx;
      wLl=(vx-xL)/(vx-vxm1);
      wLr=1-wLl;
      if(!right){
        vxp1=xR;
        wRNm1=(vxm1-xR)/(llg.vx-vx);
        wRNm2=(xR-vx)/(llg.vx-vx);
      } // !right
    }   // if left
    if(right){
      TRACE(15,"Right present");
      const LocalGeom rlg=LocalGeom(*right);
      vSfR=rlg.vSf;
      vxp1=rlg.vx;
      wRr=(xR-vx)/(vxp1-vx);
      wRl=1-wRr;
      if(!left){
        vxm1=xL;
        wL0=vxp1/(vxp1-vx);
        wL1=-vx/(rlg.vx-vx);
      } // !left
    }   // if right
  } // WeightFactors()
  void WeightFactors::show() const {
    LocalGeom::show();
    cout <<"Showing WeightFactors data..\n";    
    cout << "wLl   :"<<wLl<<"\n";
    cout << "wLr   :"<<wLr<<"\n";
    cout << "wRl   :"<<wRl<<"\n";
    cout << "wRr   :"<<wRr<<"\n";
    cout << "Special weight function:\n";
    cout << "wL0   :"<<wL0<<"\n";
    cout << "wL1   :"<<wL1<<"\n";
    cout << "wRNm1 :"<<wRNm1<<"\n";
    cout << "wRNm2 :"<<wRNm2<<"\n";
    cout << "vSfL  :"<<vSfL<<"\n";
    cout << "vSfR  :"<<vSfR<<"\n";
    cout << "SfL   :"<<SfL<<"\n";
    cout << "SfR   :"<<SfR<<"\n";
  }
} // namespace tube