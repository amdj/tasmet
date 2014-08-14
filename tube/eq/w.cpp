#include "w.h"
#include "tubevertex.h"
namespace W{
  SPOILNAMESPACE
  using segment::LocalGeom;
  W::W(){}
  void W::show() const {
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
    cout << "UsignL:"<<UsignL<<"\n";
    cout << "UsignR:"<<UsignR<<"\n";
    
    
  }
  void W::operator()(const tube::TubeVertex& v){
    const segment::Geom& geom=*v.lg.geom;
    const segment::LocalGeom& lg=v.lg;
    us i=v.i;
    wLl=wRr=wLr=wRl=0; 
    wL0=wL1=wRNm1=wRNm2=0;    	// Special boundary weight
    // functions
    UsignR=UsignL=1;

    xvim1=xvi=xvip1=0;
    vSfR=vSfL=0;
    dxm=dxp=0;

    xvi=lg.xvi;
    vSf=lg.vSf;
    if(i>0) {   
      const LocalGeom llg=geom.localGeom(i-1);
      xvim1=llg.xvi;
      dxm=xvi-xvim1;
      wLl=(lg.xvi-lg.xL)/(lg.xvi-llg.xvi);
      wLr=(lg.xL-llg.xvi)/(lg.xvi-llg.xvi);
      vSfL=llg.vSf;
    }
    if(i==0){
      const LocalGeom rlg=geom.localGeom(i+1);
      vSfL=lg.SfL;
      wL0=rlg.xvi/(rlg.xvi-lg.xvi);
      wL1=-lg.xvi/(rlg.xvi-lg.xvi);
    }
    
    if(i<v.nCells-1){
      const LocalGeom rlg=geom.localGeom(i+1);
      xvip1=rlg.xvi;
      dxp=xvip1-xvi;      
      vSfR=rlg.vSf;
      wRr=(lg.xR-lg.xvi)/(rlg.xvi-lg.xvi);
      wRl=(rlg.xvi-lg.xR)/(rlg.xvi-lg.xvi);
    }
    if(i==v.nCells-1){
      const LocalGeom llg=geom.localGeom(i-1);
      wRNm1=(llg.xvi-lg.xR)/(llg.xvi-lg.xvi);
      wRNm2=(lg.xR-lg.xvi)/(llg.xvi-lg.xvi);
      vSfR=lg.SfR;
    }    
      
  } // W(thisseg)

}	// namespace W

