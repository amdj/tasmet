#include "tubevertex.h"
#include "tube.h"

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

    TRACE(0,"TubeVertex constructor done");
  }
  TubeVertex::TubeVertex(const TubeVertex& told):TubeVertex(told.tube,told.i){
    TRACE(0,"TubeVertex::operator(),tgp");
    TRACE(-1,"Copied TubeVertex i:"<<i);
    
  }
  // void TubeVertex::operator=(const TubeVertex& rhs){
  //   TRACE(5,"Error: no copies allowed of TubeVertex");
  //   exit(1);
  // }
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
