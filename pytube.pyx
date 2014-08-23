include "pytube.pxi"
import numpy as np
cimport numpy as np


cdef extern from "models.h" namespace "":
    Solver* Fubini(us gp,us Nf,d freq,d L,d S,vd p1,int loglevel,d kappa)
    Solver* Fubini_fullenergy(us gp,us Nf,d freq,d L,d S,vd p1,int loglevel,d kappa)
    Solver* ThreeTubes(us gp,us Nf,d freq,d L,d S1,d S2,vd p1,int loglevel,d kappa)
    Solver* ConeTube(us gp,us Nf,d freq,d L,d r1,d r2,vd p1,int loglevel,d kappa)
    Solver* ThreeTubesConduction(us gp,us Nf,d freq,d L,d S1,d S2,vd p1,int loglevel,d kappa,d Tr)

    Solver* ThreeTubesEngine(us gp,us Nf,d freq,d Tr,int loglevel,d kappa)    
    Solver* ThreeTubesEngineDriven(us gp,us Nf,d freq,d Tr,vd p1,int loglevel,d kappa)    
	#Pytube compatible with paper_1 Fubini code
    
cdef class pytube:
    cdef Solver* sol
    cdef Tube* tube[5]
    cdef us ntubes
    cdef us Nf
    cdef d freq
    def __cinit__(self,us gp,us Nf,d freq,d L,d S,\
                  n.ndarray[n.float64_t,ndim=1] p1, cshape,int loglevel,d kappa,case='fubini',d S2=0,d Tr=0):
        self.sol=NULL
        self.ntubes=0
        self.Nf=Nf
        self.freq=freq
        assert(L>0)
        assert(gp>3)
        if case =='fubini':
            print('Case Fubini')
            self.sol=Fubini(gp,Nf,freq,L,S,dndtovec(p1),loglevel,kappa)
            self.tube[0]=<Tube*> self.sol[0].sys().getSeg(0)
            self.ntubes=1
        elif case =='fubini_fullenergy':
            print('Case Fubini_fullenergy')
            self.sol=Fubini_fullenergy(gp,Nf,freq,L,S,dndtovec(p1),loglevel,kappa)
            self.tube[0]=<Tube*> self.sol[0].sys().getSeg(0)
            self.ntubes=1
        elif case == 'cone':
            assert(S2>0)
            print('Case ConeTube')
            r1=np.sqrt(S/np.pi)
            r2=np.sqrt(S2/np.pi)
            self.sol=ConeTube(gp,Nf,freq,L,r1,r2,dndtovec(p1),loglevel,kappa)
            self.tube[0]=<Tube*> self.sol[0].sys().getSeg(0)
            self.ntubes=1
        elif case == 'threetubes':
            assert(S2>0)
            print('Case threetubes')            
            self.sol=ThreeTubes(gp,Nf,freq,L,S,S2,dndtovec(p1),loglevel,kappa)
            self.tube[0]=<Tube*> self.sol[0].sys().getSeg(0)
            self.tube[1]=<Tube*> self.sol[0].sys().getSeg(1)
            self.tube[2]=<Tube*> self.sol[0].sys().getSeg(2)
            self.ntubes=3
        elif case == 'threetubesenginedriven':
            assert(S2>0)
            print('Case threetubesenginedriven')            
            self.sol=ThreeTubesEngineDriven(gp,Nf,freq,Tr,dndtovec(p1),loglevel,kappa)
            self.tube[0]=<Tube*> self.sol[0].sys().getSeg(0)
            self.tube[1]=<Tube*> self.sol[0].sys().getSeg(1)
            self.tube[2]=<Tube*> self.sol[0].sys().getSeg(2)
            self.ntubes=3
        elif case == 'threetubesengine':
            assert(S2>0)
            print('Case threetubesengine')            
            self.sol=ThreeTubesEngine(gp,Nf,freq,Tr,loglevel,kappa)
            self.tube[0]=<Tube*> self.sol[0].sys().getSeg(0)
            self.tube[1]=<Tube*> self.sol[0].sys().getSeg(1)
            self.tube[2]=<Tube*> self.sol[0].sys().getSeg(2)
            self.ntubes=3
        elif case == 'threetubesconduction':
            print('Case threetubesconduction')            
            assert(S2>0)
            assert(Tr>0)
            self.sol=ThreeTubesConduction(gp,Nf,freq,L,S,S2,dndtovec(p1),loglevel,kappa,Tr)
            self.tube[0]=<Tube*> self.sol[0].sys().getSeg(0)
            self.tube[1]=<Tube*> self.sol[0].sys().getSeg(1)
            self.tube[2]=<Tube*> self.sol[0].sys().getSeg(2)
            self.ntubes=3
        else:
            print('Warning: no valid case selected! Tried was: %s' %case)

    def show(self,showvertex):
        if showvertex is False:
            self.sol[0].sys().show(False)
        else:
            self.sol[0].sys().show(True)
    def __dealloc__(self):
        if self.sol!=NULL:
            del self.sol
    cpdef getNf(self):
        return self.Nf
    cpdef getFreq(self):
        return self.freq
    cpdef solve(self,maxiter=10):
        self.sol[0].solve()    
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
    cpdef setRes(self,n.ndarray[n.float64_t,ndim=1] res):
        self.sol[0].sys().setRes(dndtovec(res))

    cpdef getErrorEq(self,eqnr,freqnr,i=0):
        assert(eqnr<5)
        assert(freqnr<2*self.Nf+1)
        assert(i<self.ntubes)
        return dvectond(self.tube[i].getErrorAt(eqnr,freqnr))
    
    cpdef getResVar(self,_type,freqnr,i=0):
        assert(i<self.ntubes)
        if _type=='pres':
            return dvectond(self.tube[i].getResAt(3,freqnr))
        elif _type=='rho':
            return dvectond(self.tube[i].getResAt(0,freqnr))
        elif _type=='temp':
            return dvectond(self.tube[i].getResAt(2,freqnr))
        elif _type=='volu':
            return dvectond(self.tube[i].getResAt(1,freqnr))
        elif _type=='stemp':
            return dvectond(self.tube[i].getResAt(4,freqnr))
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

