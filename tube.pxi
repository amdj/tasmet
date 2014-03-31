include "armapy.pxi"

cdef extern from "logger.h" namespace "":
    void initlog(us)

cdef extern from "globalconf.h" namespace "tasystem":
    cdef cppclass Globalconf:
        Globalconf(us Nf,d freq)
        Globalconf(us Nf,d freq,string Gas,d T0,d p0,d Mass)

cdef extern from "tube/geom.h" namespace "tube":
    cdef cppclass Geom:
        Geom(us gp,d L,d S,d phi,d rh,string cshape) except +
        vd vx


cdef extern from "var/var.h" namespace "variable":
    cdef cppclass varoperations:
        varoperations(us,d)
    cdef cppclass var:
        var(varoperations&,double)   #Initialize with constant value
        void set(double,us) #Set frequency us to double

        
cdef extern from "tube/tube.h" namespace "tube":
    cdef cppclass Tube:
        Tube(Globalconf&,Geom g)
        void Init(d T0,d p0)
        us Ncells
        Geom geom
        void setLeftbc(Vertex* v)    
        varoperations vop    
        vd Error()
        vd GetRes()
        void DoIter()
        vd GetResAt(us varnr,us freqnr)
cdef extern from "tube/vertex.h":
    cdef cppclass Vertex:
        pass
    cdef cppclass TubeVertex(Vertex):
        pass

cdef extern from "tube/bcvertex.h" namespace "tube":
    cdef cppclass LeftPressure(TubeVertex):
        LeftPressure(Tube,var pres,d T0)
            
    

