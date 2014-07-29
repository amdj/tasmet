include "pytube.pxi"
import numpy as np
cimport numpy as np


cdef extern from "fubini.h" namespace "":
    Solver* Fubini(us gp,us Nf,d freq,d L,d S,vd p1,int loglevel,d kappa)
cdef extern from "threetubes.h" namespace "":
    Solver* ThreeTubes(us gp,us Nf,d freq,d L,d S1,d S2,vd p1,int loglevel,d kappa)
cdef extern from "conetube.h" namespace "":
    Solver* ConeTube(us gp,us Nf,d freq,d L,d r1,d r2,vd p1,int loglevel,d kappa)

    
#Pytube compatible with paper_1 Fubini code
    
cdef class pytube:
    cdef Solver* sol
    cdef Tube* tube[5]
    cdef us ntubes
    def __cinit__(self,us gp,us Nf,d freq,d L,d S,d T0,d p0,\
                  n.ndarray[n.float64_t,ndim=1] p1, cshape,int loglevel,d kappa,case='fubini',d S2=0):
        self.sol=NULL
        self.ntubes=0
        assert(L>0)
        assert(gp>3)
        if case =='fubini':
            print('Case Fubini')
            self.sol=Fubini(gp,Nf,freq,L,S,dndtovec(p1),loglevel,kappa)
            self.tube[0]=<Tube*> self.sol[0].sys[0].getSeg(0)
            self.ntubes=1
        elif case == 'cone':
            assert(S2>0)
            r1=np.sqrt(S/np.pi)
            r2=np.sqrt(S2/np.pi)
            self.sol=ConeTube(gp,Nf,freq,L,r1,r2,dndtovec(p1),loglevel,kappa)
            self.tube[0]=<Tube*> self.sol[0].sys[0].getSeg(0)
            self.ntubes=1
        elif case == 'threetubes':
            assert(S2>0)
            
            self.sol=ThreeTubes(gp,Nf,freq,L,S,S2,dndtovec(p1),loglevel,kappa)
            self.tube[0]=<Tube*> self.sol[0].sys[0].getSeg(0)
            self.tube[1]=<Tube*> self.sol[0].sys[0].getSeg(1)
            self.tube[2]=<Tube*> self.sol[0].sys[0].getSeg(2)
            self.ntubes=3
        else:
            print('Warning: no valid case selected! Tried was: %s' %case)

    def show(self):
        self.sol[0].sys[0].show()    
    def __dealloc__(self):
        if self.sol!=NULL:
            del self.sol
    cpdef getNf(self):
        return self.Nf
    cpdef solve(self):
        self.sol[0].solve()    
    cpdef getgp(self):
        return self.gp
    cpdef getx(self,i=0):
        assert(i<self.ntubes)
        return dvectond(self.tube[i].geom.vx)
    cpdef getSf(self,i=0):
        assert(i<self.ntubes)
        return dvectond(self.tube[i].geom.vSf)

    cpdef Error(self):
        return dvectond(self.sol[0].sys[0].Error())
    cpdef GetRes(self):
        return dvectond(self.sol[0].sys[0].GetRes())
    cpdef DoIter(self,d relaxfac):
        self.sol[0].DoIter(relaxfac)
    cpdef SetRes(self,n.ndarray[n.float64_t,ndim=1] res):
        self.sol[0].sys[0].SetRes(dndtovec(res))
    cpdef GetResVar(self,_type,freqnr,i=0):
        assert(i<self.ntubes)
        if _type=='pres':
            return dvectond(self.tube[i].GetResAt(3,freqnr))
        elif _type=='rho':
            return dvectond(self.tube[i].GetResAt(0,freqnr))
        elif _type=='temp':
            return dvectond(self.tube[i].GetResAt(2,freqnr))
        elif _type=='volu':
            return dvectond(self.tube[i].GetResAt(1,freqnr))
        else:
            return None
    # def getResult(self):
    #     cdef vector[double] res=self.thisl.getRes()
    #     return ToArray(res)
    # def getError(self):
    #     cdef vector[double] er=self.thisl.getErr()
    #     return ToArray(er)
    # def getpResult(self,int j):
    #     cdef vector[double] res=self.thisl.p.getRes(j)
    #     return ToArray(res)
    # def getUResult(self,int j):
    #     cdef vector[double] res=self.thisl.U.getRes(j)
    #     return ToArray(res)
    # def getrhoResult(self,int j):
    #     cdef vector[double] res=self.thisl.rho.getRes(j)
    #     return ToArray(res)
    # def setResult(self,n.ndarray[n.float64_t,ndim=1] res):
    #     self.thisl.setResult(ToVec(res))
    # def setpResult(self,n.ndarray[n.float64_t,ndim=1] res,int f):
    #     self.thisl.p.setResult(ToVec(res),f)
    # def setUResult(self,n.ndarray[n.float64_t,ndim=1] res,int f):
    #     self.thisl.U.setResult(ToVec(res),f)
    # def Objfun(self,n.ndarray[n.float64_t,ndim=1] guess):
    #     self.setResult(guess)
    #     return self.getError()
        # del self.tube1

