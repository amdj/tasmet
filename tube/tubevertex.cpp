#include "tubevertex.h"
#include "tube.h"
#include "pressurebc.h"
namespace tube{

  TubeVertex::TubeVertex(const Tube& tube1,us i):Vertex(tube1,i),tube(tube1),c(tube,*this),m(tube,*this),e(tube,*this),s(tube,*this),se(tube,*this),is(tube,*this)
  {
    TRACE(0,"TubeVertex contructor");
    
    eq[0]=&this->c;			// Continuity is first
    eq[1]=&this->m;
    eq[2]=&is; 			// Changed to isentropic
    // eq[2]=&e; 			// Full energy
    eq[3]=&s;
    eq[4]=&se;
  }

  TubeVertex::TubeVertex(const TubeVertex& told):TubeVertex(told.tube,told.i){
    TRACE(0,"TubeVertex::operator(),tgp");
    TRACE(-1,"Copied TubeVertex i:"<<i);
  }
  void TubeVertex::updateW(){
    TRACE(1,"TubeVertex::updateW()");
    Vertex::updateW();
    const us& Ncells=tube.Ncells;
    const Geom& geom=tube.geom;
    
    const vd& vx=tube.geom.vx;
    const d& vxi=vx(i);
    d vxip1=0;
    d vxim1=0;
    // Initialize distances to next node to zero
    dxm=dxp=0;

    // Left and right cross-sectional area
    SfL=geom.Sf(i);
    SfR=geom.Sf(i+1);
    // Geometric parameters
    vSf=geom.vSf(i);
    vSs=geom.vSs(i);
    vVf=geom.vVf(i);
    vVs=geom.vVs(i);

    xR=tube.geom.x(i+1);		// Position of right cell wall
    xL=tube.geom.x(i);			// Position of left cell wall

    int UsignL=1;    
    int UsignR=1;    

    // Initialize weight functions to zero
    wLl=0; wLr=0; wRr=0; wRl=0;
    wL0=wL1=wRNm1=wRNm2=0;	// Put these weight functions to zero

    // ****************************** Initalization of vxipm and dxpm
    if(i>0){
      vxim1=vx(i-1);
      dxm=vxi-vxim1;
      // Left weight functions
      wLl=(vxi-xL)/(vxi-vxim1);
      wLr=(xL-vxim1)/(vxi-vxim1);
    }
    if(i<Ncells-1){
      vxip1=vx(i+1);
      dxp=vxip1-vxi;
      // Right weight functions
      wRr=(xR-vxi)/(vxip1-vxi);
      wRl=(vxip1-xR)/(vxip1-vxi);
    }
    // ****************************** End initialization
    // ****************************** Initialization special weight functions
    if(i==Ncells-1){
      wRNm1=(vxim1-xR)/(vxim1-vxi);
      wRNm2=(xR-vxi)/(vxim1-vxi);
    }
    if(i==0){
      wL0=vxip1/(vxip1-vxi);
      wL1=-vxi/(vxip1-vxi);
    }
    // ****************************** End Special weight functions

    if(i==0 && left!=NULL){
      if(left->seg.gettype().compare("Tube")==0){ // Its a Tube
	if(*left->seg.right==seg){
	  TRACE(5,"Connected current head to left segment's tail");
	  const us& LeftNcells=left->seg.Ncells;
	  d L=left->seg.geom.L;
	  vxim1=left->seg.geom.vx(LeftNcells-1)-L;
	  TRACE(6,"vxi:"<<vxi);
	  TRACE(6,"vxim1:"<<vxim1);
	}
      	else{
	  TRACE(5,"Connected current head to left segment's head");
	  TRACE(6,"vxi:"<<vxi);
	  TRACE(6,"vxim1:"<<vxim1);
	  vxim1=-left->seg.geom.vx(0);
	  wLl=(vxi-xL)/(vxi-vxim1);
	  UsignL=-1;
      	}
	wLl=(vxi-xL)/(vxi-vxim1);	
	wLr=(xL-vxim1)/(vxi-vxim1);
      }
      else{			// Its not a Tube
	TRACE(20,"Error, this kind of coupling not yet implemented. Exiting...");
	TRACE(20,"Seg on left side: "<< right->seg.gettype());
	exit(1);
      }
    } // i==0 and left!=NULL
    
    else if((i==Ncells-1) && right!=NULL) {
      if(right->seg.gettype().compare("Tube")==0){ // Its a Tube
	if(*right->seg.left==seg){
	  TRACE(5,"Connected current tail to right segment's head");
	  d L=seg.geom.L;
	  vxip1=right->seg.geom.vx(0)+seg.geom.L;
	  TRACE(6,"vxi:"<<vxi);
	  TRACE(6,"vxip1:"<<vxip1);
	}
      	else{
	  TRACE(5,"Connected current tail to right segment's tail");
	  const us& RightNcells=right->seg.Ncells;
	  const d& RightL=right->seg.geom.L;
	  vxip1=RightL-right->seg.geom.vx(Ncells)+seg.geom.L;
	  TRACE(6,"vxi:"<<vxi);
	  TRACE(6,"vxip1:"<<vxip1);
	  UsignR=-1;
      	}
	wRr=(xR-vxi)/(vxip1-vxi);
	wRl=(vxip1-xR)/(vxip1-vxi);
      }
      else{			// Its not a Tube
	TRACE(20,"Error, this kind of coupling not yet implemented. Exiting...");
	TRACE(20,"Seg on left side: "<< right->seg.gettype());
	exit(1);
      }
    } // i==Ncells-1 and right!=NULL

    c.Wddt=vVf;
    m.Wddt=vVf/vSf;
    e.Wddt=vVf;

    if( (i>0 && i<Ncells-1) || (i==0 && left!=NULL) || (i==Ncells-1 && right!=NULL) ) {
      c.Wim1=-UsignL*wLl;
      c.Wi=UsignR*wRl-UsignL*wLr;
      c.Wip1=UsignR*wRr;

      m.Wuim1=-wLl/SfL;
      m.Wui=wRl/SfR-wLr/SfL;
      m.Wuip1=wRr/SfR;

      m.Wpim1=-SfL*wLl;
      m.Wpi  = SfR*wRl-SfL*wLr;
      m.Wpip1= SfR*wRr;

      e.Wgim1=-wLl;
      e.Wgi=wRl-wLr;
      e.Wgip1=wRr;

      e.Wjim1=wLl;
      e.Wji=wLr-wRl;
      e.Wjip1=-wRr;

      e.Wc1=-SfL/dxm;
      e.Wc2=SfL/dxm;
      e.Wc3=SfR/dxp;
      e.Wc4=-SfR/dxp;
    }
    else if(i==0){
      c.Wim1=0;
      c.Wi=wRl;
      c.Wip1=wRr;
      
      m.Wuim1=0;
      m.Wui=wRl/SfR;
      m.Wuip1=wRr/SfR;
      
      m.Wpim1=0;
      m.Wpi=SfR*wRl-SfL*wL0;
      m.Wpip1=SfR*wRr-SfL*wL1;

      e.Wgim1=0;
      e.Wgi=wRl;
      e.Wgip1=wRr;

      e.Wjim1=0;
      e.Wji=wL0-wRl;
      e.Wjip1=wL1-wRr;

      e.Wc1=0;
      e.Wc2=0;
      e.Wc3=SfR/dxp;
      e.Wc4=-SfR/dxp;
    }
    else if(i==Ncells-1){
      c.Wi=-wLr;
      c.Wim1=-wLl;
      c.Wip1=0;

      m.Wuim1=-wLl/SfL;
      m.Wui=-wLr/SfL;
      m.Wuip1=0;
      
      m.Wpim1=-SfL*wLl+SfR*wRNm2;
      m.Wpi=-SfL*wLr+SfR*wRNm1;
      m.Wpip1=0;

      e.Wgim1=-wLl;
      e.Wgi=-wLr;
      e.Wgip1=0;

      e.Wjim1=wLl-wRNm2;
      e.Wji=wLr-wRNm1;
      e.Wjip1=0;

      e.Wc1=-SfL/dxm;
      e.Wc2=SfL/dxm;
      e.Wc3=0;
      e.Wc4=0;
    } 

    // Contribution from changing cross-sectional area
    m.Wpi+=SfL-SfR;

  }
  vd  TubeVertex::csource() const {
    TRACE(0,"TubeVertex::csource()");
    return zeros(Ns);}
  vd  TubeVertex::msource() const {
    TRACE(0,"TubeVertex::msource()");
    return zeros(Ns);}
  vd  TubeVertex::esource() const {
    TRACE(0,"TubeVertex::esource()");
    return zeros(Ns);}    
    
  TubeVertex::~TubeVertex(){
    TRACE(-5,"TubeVertex destructor");
  }

 



} // namespace tube
