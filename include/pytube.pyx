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
        return self.sol[0].sys().gc.Nf()
    cpdef getL(self,nr):
        assert(nr<self.ntubes)
        return self.tube[nr].geom().L()
    cpdef nTubes(self):
        return self.ntubes
    cpdef setNf(self,Nf):
        self.sol.sys().setNf(Nf)
    cpdef getFreq(self):
        return self.sol[0].sys().gc.getfreq()
    cpdef setFreq(self,d freq):
        self.sol[0].sys().gc.setfreq(freq)
    cpdef solve(self,us maxiter=4000,d mindampfac=0.1,d maxdampfac=1.0,d funtol=1e-9,d reltol=1e-6,wait=True):
        self.sol[0].solve(maxiter,funtol,reltol,mindampfac,maxdampfac,wait)
    cpdef stop(self):
        self.sol[0].stop()    
    cpdef cls(self):
        clearConsole()
    cpdef getx(self,i=0):
        assert(i<self.ntubes)
        return dvectond(self.tube[i].geom().vx_vec())
    cpdef getS(self,i=0):
        assert(i<self.ntubes)
        return dvectond(self.tube[i].geom().vS_vec())
    cpdef getphi(self,i=0):
        assert(i<self.ntubes)
        return dvectond(self.tube[i].geom().vphi_vec())
    cpdef getrh(self,i=0):
        assert(i<self.ntubes)
        return dvectond(self.tube[i].geom().vrh_vec())
    cpdef getSf(self,i=0):
        assert(i<self.ntubes)
        return dvectond(self.tube[i].geom().vSf_vec())
    cpdef error(self):
        return eigentond(self.sol[0].sys().error())
    cpdef getRes(self):
        return eigentond(self.sol[0].sys().getRes())
    cpdef doIter(self,d relaxfac=1.0):
        self.sol[0].doIter(relaxfac)
    cpdef setRes(self,res):
        self.sol[0].sys().setRes(dndtovec(res))
    cpdef setResPt(self, pytubeBase other):
        self.sol[0].sys().setRes(other.sol[0].sys())
    cpdef resetHarmonics(self):
        self.sol[0].sys().resetHarmonics()
    
    cpdef getErrorEq(self,eqnr,freqnr,tubenr=0):
        assert(eqnr<5)
        assert(freqnr<2*self.getNf()+1)
        assert(tubenr<self.ntubes)
        return dvectond(self.tube[tubenr].getErrorAt(eqnr,freqnr))

    cpdef getHtot(self,tubenr):    
        i=tubenr
        assert(i<self.ntubes)
        return dvectond(self.tube[i].Htot())
    cpdef dragCoef(self,freqnr,tubenr=0):
        return dvectond(self.tube[tubenr].dragCoefVec(freqnr))
    cpdef getResVar(self,_type,us freqnr,us tubenr=0):
    # Unfortunately, this is hard-coded
        i=tubenr
        assert(i<self.ntubes)
        assert(freqnr<(2*self.getNf()+1))
        j=0
        if _type=='rho':
            j=0
        elif _type=='volu':
            j=1
        elif _type=='temp':
            j=2
        elif _type=='pres':
            j=3
        elif _type=='stemp':
            j=4
        else:
            return None
        return dvectond(self.tube[i].getResAt(j,freqnr))
    cpdef showJac(self):
        self.sol[0].sys().showJac()

        
