include "armapy.pxi"

cdef extern from "logger.h" namespace "":
    void initlog(us)

cdef extern from "globalconf.h" namespace "tasystem":
    cdef cppclass Globalconf:
        Globalconf(us Nf,d freq,string Gas,d T0,d p0,d Mach,d S,d dx,d Mass,d kappa)

cdef extern from "tube/geom.h" namespace "tube":
    cdef cppclass Geom:
        Geom(us gp,d L,d S,d phi,d rh,string cshape) except +
        vd vx

cdef extern from "var/var.h" namespace "variable":
    cdef cppclass var:
        var(Globalconf&,double)   #Initialize with constant value
        var(Globalconf&)   #Initialize with zeros
        void set(double,us) #Set frequency us to double

        
cdef extern from "tube/tube.h" namespace "tube":
    cdef cppclass Tube:
        Tube(Globalconf&,Geom g)
        void Init(d T0,d p0)
        us Ncells
        Geom geom
        Globalconf& gc
        void setLeftbc(Vertex* v)
        void setRightbc(Vertex* v)                
        vd Error()
        vd GetRes()
        void SetRes(vd res)    
        void DoIter(d)
        vd GetResAt(us varnr,us freqnr)
cdef extern from "tube/vertex.h":
    cdef cppclass Vertex:
        pass
    cdef cppclass TubeVertex(Vertex):
        pass

cdef extern from "tube/bcvertex.h" namespace "tube":
    cdef cppclass LeftPressure(TubeVertex):
        LeftPressure(Tube,var pres)
    cdef cppclass RightImpedance(TubeVertex):        
        RightImpedance(Tube,vd Z)    















