#include "w.h"
#include "tubevertex.h"
namespace W{
  SPOILNAMESPACE
  using segment::LocalGeom;
  W::W(){}
  void W::operator()(const tube::TubeVertex& v){
    const segment::Geom& geom=*v.lg.geom;
    const segment::LocalGeom& lg=v.lg;
    us i=v.i;
    wLl=wRr=wLr=wRl=0; 
    wL0=wL1=wRNm1=wRNm2=0;    	// Special boundary weight
    // functions
    UsignR=UsignL=1;

    vxim1=vxi=vxip1=0;
    vSfR=vSfL=0;
    dxm=dxp=0;

    vxi=lg.vxi;
    if(i>0) {   
      const LocalGeom llg=geom.localGeom(i-1);
      vxim1=llg.vxi;
      dxm=vxi-vxim1;
      wLl=(lg.vxi-lg.xL)/(lg.vxi-llg.vxi);
      wLr=(lg.xL-llg.vxi)/(lg.vxi-llg.vxi);
      vSfL=llg.vSf;
      dxm=vxi-vxim1;
    }
    if(i==0){
      const LocalGeom rlg=geom.localGeom(i+1);
      wL0=rlg.vxi/(rlg.vxi-lg.vxi);
      wL1=-rlg.vxi/(rlg.vxi-lg.vxi);
    }
    
    if(i<v.nCells-1){
      const LocalGeom rlg=geom.localGeom(i+1);
      vxip1=rlg.vxi;
      dxp=vxip1-vxi;      
      wRr=(lg.xR-lg.vxi)/(rlg.vxi-lg.vxi);
      wRl=(rlg.vxi-lg.xR)/(rlg.vxi-lg.vxi);
      vSfR=rlg.vSf;
      dxp=vxip1-vxi;
    }
    if(i==v.nCells-1){
      const LocalGeom llg=geom.localGeom(i-1);
      wRNm1=(llg.vxi-lg.xR)/(llg.vxi-lg.vxi);
      wRNm2=(lg.xR-lg.vxi)/(llg.vxi-lg.vxi);
    }    
      
  } // W(thisseg)

}	// namespace W

