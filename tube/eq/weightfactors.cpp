#include "weightfactors.h"

namespace tube{



  WeightFactors::WeightFactors(const TubeVertex& v):
    LocalGeom(v.getTube().geom(),v.geti())
  {
    TRACE(15,"WeightFactors::WeightFactors()");
    const TubeVertex* left=v.left();
    const TubeVertex* right=v.right();

    if(left) {   
      const LocalGeom& llg=left->localGeom();
      vxm1=llg.vx;
      dxm=vx-vxm1;
      wLl=(vx-xL)/(vx-llg.vx);
      wLr=(xL-llg.vx)/(vx-llg.vx);
      vSfL=llg.vSf;
    }
    else{
      const LocalGeom& rlg=right->localGeom();
      vSfL=SfL;
      wL0=rlg.vx/(rlg.vx-vx);
      wL1=-lg.vx/(rlg.vx-vx);
    }
    if(right){
      const LocalGeom& rlg=right->lg;
      vxp1=rlg.vx;
      dxp=vxp1-vx;      
      vSfR=rlg.vSf;
      wRr=(xR-vx)/(rlg.vx-vx);
      wRl=(rlg.vx-xR)/(rlg.vx-vx);
    }
    else{
      const LocalGeom& llg=left->lg;
      wRNm1=(llg.vx-xR)/(llg.vx-vx);
      wRNm2=(lg.xR-vx)/(llg.vx-vx);
      vSfR=SfR;
    }    

  } // WeightFactors()
  void WeightFactors::show() const {
    LocalGeom::show();
    cout <<"Showing WeightFactors data..\n";    
    cout << "wLl   :"<<wLl<<"\n";
    cout << "wLr   :"<<wLr<<"\n";
    cout << "wRl   :"<<wRl<<"\n";
    cout << "wRr   :"<<wRr<<"\n";
    cout << "wL0   :"<<wL0<<"\n";
    cout << "wL1   :"<<wL1<<"\n";
    cout << "wRNm1 :"<<wRNm1<<"\n";
    cout << "wRNm2 :"<<wRNm2<<"\n";
    cout << "vSfL  :"<<vSfL<<"\n";
    cout << "vSfR  :"<<vSfR<<"\n";
    cout << "dxm   :"<<dxm<<"\n";
    cout << "dxp   :"<<dxp<<"\n";
  }
} // namespace tube
