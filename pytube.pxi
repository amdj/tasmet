include "armapy.pxi"

cdef extern from "logger.h" namespace "":
    void initlog(us)

cdef extern from "globalconf.h" namespace "tasystem":
    cdef cppclass Globalconf:
        Globalconf(us Nf,d freq,string Gas,d T0,d p0,d Mach,d S,d dx,d Mass,d kappa)

# cdef extern from <vector.h> namespace "std":
    # cdef cppclass vector[T]:
        # T operator[](int)
        
cdef extern from "seg.h" namespace "segment":
    cdef cppclass Seg:
        Seg(Globalconf& gc)

        
cdef extern from "system.h" namespace "tasystem":
    cdef cppclass TAsystem:
        TAsystem(Globalconf& gc)
        void addseg(Seg& seg)
        vd Error()
        void show()
        vd GetRes()
        void SetRes(vd)
        Seg* getSeg(us i)
        
cdef extern from "solver.h" namespace "tasystem":
    cdef cppclass Solver:
        Solver(TAsystem& sys)
        void DoIter(d dampfac)
        TAsystem* sys
        void Init()
    

cdef extern from "geom.h" namespace "tube":
    cdef cppclass Geom:
        Geom(us gp,d L,d S,d phi,d rh,string cshape) except +
        vd vx

cdef extern from "var.h" namespace "variable":
    cdef cppclass var:
        var(Globalconf&,double)   #Initialize with constant value
        var(Globalconf&)   #Initialize with zeros
        void set(double,us) #Set frequency us to double
        void set(vd) #Set frequency us to double


cdef extern from "tube.h" namespace "tube":
    cdef cppclass Tube:
        Tube(Geom g)
        Seg(Geom g)    
        void Init(d T0,d p0)
        us Ncells
        Geom geom
        Globalconf& gc
        void setLeftbc(Vertex* v)
        void setRightbc(Vertex* v)                
        vd Error()
        vd GetRes()
        void SetRes(vd res)    
        vd GetResAt(us varnr,us freqnr)
cdef extern from "vertex.h":
    cdef cppclass Vertex:
        pass
    cdef cppclass TubeVertex(Vertex):
        pass
cdef extern from "bcvertex.h" namespace "tube":
    cdef cppclass TubeVertex:
         TubeVertex(Tube&,us i)
cdef extern from "pressurebc.h" namespace "tube":
    cdef cppclass LeftPressure(TubeVertex):
        LeftPressure(Tube,var pres)
cdef extern from "impedancebc.h" namespace "tube":
    cdef cppclass RightImpedance(TubeVertex):        
        RightImpedance(Tube,vd Z)    















