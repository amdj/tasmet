include "pytube.pxi"
import numpy as np
cimport numpy as np


cdef class pytubeBase:
    cdef Solver* sol
    cdef Tube* tube[20]
    cdef us ntubes
    def show(self,showvertex=False):
        if showvertex is False:
            self.sol[0].sys().show(0)
        else:
            self.sol[0].sys().show(1)
    def __dealloc__(self):
        if self.sol!=NULL:
            del self.sol
    cpdef getNf(self):
        return self.sol[0].sys().gc.Nf
    cpdef getL(self,nr):
        assert(nr<self.ntubes)
        return self.tube[nr].geom.L
    cpdef nTubes(self):
        return self.ntubes
    cpdef getFreq(self):
        return self.sol[0].sys().gc.getfreq()
    cpdef solve(self,us maxiter=100,d dampfac=1.0,d funtol=1e-9,d reltol=1e-6):
        self.sol[0].solve(maxiter,funtol,reltol,dampfac)    
    cpdef getgp(self):
        return self.gp
    cpdef getx(self,i=0):
        assert(i<self.ntubes)
        return dvectond(self.tube[i].geom.xv)
    cpdef getSf(self,i=0):
        assert(i<self.ntubes)
        return dvectond(self.tube[i].geom.vSf)
    cpdef error(self):
        return eigentond(self.sol[0].sys().error())
    cpdef getRes(self):
        return eigentond(self.sol[0].sys().getRes())
    cpdef doIter(self,d relaxfac=1.0):
        self.sol[0].doIter(relaxfac)
    cpdef setRes(self,res):
        self.sol[0].sys().setRes(dndtovec(res))

    cpdef setResPt(self, pytubeBase other):
        self.sol.sys().setRes(other.sol[0].sys())

    
    cpdef getErrorEq(self,eqnr,freqnr,tubenr=0):
        assert(eqnr<5)
        assert(freqnr<2*self.getNf()+1)
        assert(tubenr<self.ntubes)
        return dvectond(self.tube[tubenr].getErrorAt(eqnr,freqnr))
    cpdef getResVar(self,_type,freqnr,i=0):
        assert(i<self.ntubes)
        if _type=='pres':
            return dvectond(self.tube[i].getResAt(2,freqnr))
        elif _type=='rho':
            return dvectond(self.tube[i].getResAt(0,freqnr))
        elif _type=='temp':
            return dvectond(self.tube[i].getResAt(3,freqnr))
        elif _type=='volu':
            return dvectond(self.tube[i].getResAt(1,freqnr))
        elif _type=='stemp':
            return dvectond(self.tube[i].getResAt(4,freqnr))
        else:
            return None
    cpdef showJac(self):
        self.sol[0].sys().showJac()

        
cdef extern from "models.h" namespace "":
    Solver* Fubini(us gp,us Nf,d freq,d L,vd p1,int loglevel,d kappa,int options)    
    Solver* ThreeTubes(us gp,us Nf,d freq,d L,d R1,d R2,vd p1,int loglevel,d kappa,d Tr,int options)
    Solver* Atchley_Engine(us gp,us Nf,d freq,d Tr,int loglevel,d kappa,vd p1,d p0,int options) 
    Solver* SimpleTube(us gp,us Nf,d freq,d L,d r,d Tl,d Tr,vd p1,int loglevel,d kappa,int options)

cdef class simpletube(pytubeBase):        
    def __cinit__(self,us gp,us Nf,d freq,d L,d r,d Tl,d Tr
                  ,n.ndarray[n.float64_t,ndim=1] p1,\
                  int loglevel,d kappa,int options):
        self.sol=SimpleTube(gp,Nf,freq,L,r,Tl,Tr,dndtovec(p1),loglevel,kappa,options)
        self.tube[0]=<Tube*> self.sol[0].sys().getSeg(0)
        self.ntubes=1

    
cdef class fubini(pytubeBase):        
    def __cinit__(self,us gp,us Nf,d freq,d L\
                  ,n.ndarray[n.float64_t,ndim=1] p1,int loglevel,d kappa,int options):
        self.sol=Fubini(gp,Nf,freq,L,dndtovec(p1),loglevel,kappa,options)
        self.tube[0]=<Tube*> self.sol[0].sys().getSeg(0)
        self.ntubes=1

cdef class threetubes(pytubeBase):        
    def __cinit__(self,us gp,us Nf,d freq,d L,d R1,d R2\
                  ,n.ndarray[n.float64_t,ndim=1] p1,int loglevel,d kappa,d Tr,int options):
        self.sol=ThreeTubes(gp,Nf,freq,L,R1,R2,dndtovec(p1),loglevel,kappa,Tr,options)
        self.tube[0]=<Tube*> self.sol[0].sys().getSeg(0)
        self.tube[1]=<Tube*> self.sol[0].sys().getSeg(1)
        self.tube[2]=<Tube*> self.sol[0].sys().getSeg(2)
        self.ntubes=3

cdef class atchley(pytubeBase):
    def __cinit__(self,gp,Nf,freq,p1,Tr\
                  ,int loglevel,d kappa,d p0,options):
        self.sol=Atchley_Engine(gp,Nf,freq,Tr,\
                                loglevel,kappa,dndtovec(p1),p0,options)
        self.tube[0]=<Tube*> self.sol[0].sys().getSeg(0)
        self.tube[1]=<Tube*> self.sol[0].sys().getSeg(1)
        self.tube[2]=<Tube*> self.sol[0].sys().getSeg(2)
        self.tube[3]=<Tube*> self.sol[0].sys().getSeg(3)
        self.tube[4]=<Tube*> self.sol[0].sys().getSeg(4)
        self.ntubes=5

