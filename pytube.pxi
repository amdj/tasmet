include "armapy.pxi"

cdef extern from "logger.h" namespace "":
    void initlog(us)

cdef extern from "globalconf.h" namespace "tasystem":
    cdef cppclass Globalconf:
        Globalconf(us Nf,d freq,string Gas,d T0,d p0,d Mach,d S,d dx,d Mass,d kappa)

cdef extern from "seg.h" namespace "segment":
    cdef cppclass Seg:
        Seg(Globalconf& gc)

        
cdef extern from "system.h" namespace "tasystem":
    cdef cppclass taSystem:
        taSystem(Globalconf& gc)
        void addSeg(Seg& seg)
        vd error()
        void show(bool)
        vd getRes()
        void setRes(vd)
        Seg* getSeg(us i)
        
cdef extern from "solver.h" namespace "tasystem":
    cdef cppclass Solver:
        Solver(taSystem& sys)
        void doIter(d dampfac)
        taSystem& sys()
        void init()
        void solve(us maxiter)
        void solve()    

cdef extern from "geom.h" namespace "tube":
    cdef cppclass Geom:
        Geom(us gp,d L,d S,d phi,d rh,string cshape) except +
        vd xv
        vd vSf
        vd vSs
        vd x

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
        void init(d T0,d p0)
        Geom geom
        Globalconf& gc
        vd error()
        vd getRes()
        void setRes(vd res)    
        vd getResAt(us varnr,us freqnr)
        vd getErrorAt(us eqnr,us freqnr)            
cdef extern from "vertex.h":
    cdef cppclass Vertex:
        pass
cdef extern from "tubevertex.h":
    cdef cppclass TubeVertex(Vertex):
        pass















