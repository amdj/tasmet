include "armapy.pxi"
from libcpp cimport bool

cdef extern from "consolecolors.h" namespace "":
    void clearConsole()
cdef extern from "tracer.h" namespace "":
    void initlog(us)

cdef extern from "globalconf.h" namespace "tasystem":
    cdef cppclass Globalconf:
        Globalconf(us Nf,d freq,string Gas,d T0,d p0,d Mach,d S,d dx,d Mass,d kappa)
        d getfreq()
        void setfreq(d)        
        d getomg()
        void setNf(us)
        us Nf()

cdef extern from "seg.h" namespace "segment":
    cdef cppclass Seg:
        Seg(Globalconf& gc)

        
cdef extern from "tasystem.h" namespace "tasystem":
    cdef cppclass TaSystem:
        TaSystem(Globalconf& gc)
        Globalconf gc
        void addSeg(Seg& seg)
        evd error()
        void show(us)
        evd getRes()
        void setRes(const TaSystem&)
        void setRes(vd)
        Seg* getSeg(us i)
        void showJac()
        void resetHarmonics()
        void setNf(us)
        
cdef extern from "solver.h" namespace "tasystem":
    cdef cppclass Solver:
        Solver(TaSystem& sys)
        void doIter(d dampfac)
        TaSystem& sys()
        void init()
        void solve(us maxiter,d funtol,d reltol,d mindampfac,d maxdampfac,bool wait)
        void stop()


cdef extern from "geom.h" namespace "tube":
    cdef cppclass Geom:
        Geom(us gp,d L,d S,d phi,d rh,string cshape) except +
        vd vSf_vec()
        vd vrh_vec()
        vd vS_vec()
        vd vphi_vec()        
        vd x()
        vd vx_vec()
        d L()

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
        Geom geom()
        Globalconf& gc
        vd error()
        vd getRes()
        void setRes(vd res)    
        vd getResAt(us varnr,us freqnr)
        vd getErrorAt(us eqnr,us freqnr)
        vd Htot()
        vd dragCoefVec(us freqnr)
        
cdef extern from "vertex.h":
    cdef cppclass Vertex:
        pass
cdef extern from "tubevertex.h":
    cdef cppclass TubeVertex(Vertex):
        pass















