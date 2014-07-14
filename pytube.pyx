include "pytube.pxi"

cdef extern from "fubini.h" namespace "":
    Solver* Fubini(us gp,us Nf,d freq,d L,d S,vd p1,int loglevel,d kappa)

cdef class pytube:
    cdef Solver* sol
    cdef Tube* tube
    def __cinit__(self,us gp,us Nf,d freq,d L,d S,d T0,d p0,n.ndarray[n.float64_t,ndim=1] p1, cshape,int loglevel,d kappa):
        # self.thisl=new isentropictube(gp,Nf,1.,1.,freq,Up)
        # print "New tube initialized"
        # self.Nf=Nf
        # self.freq=freq
        # self.gp=gp
        self.sol=Fubini(gp,Nf,freq,L,S,dndtovec(p1),loglevel,kappa)
        self.tube=<Tube*> self.sol[0].sys[0].getSeg(0)

    def __dealloc__(self):
        del self.sol
    def init(self):
        self.sol[0].Init()    
    cpdef getNf(self):
        return self.Nf
    cpdef getgp(self):
        return self.gp
    cpdef getx(self):
        return dvectond(self.tube[0].geom.vx)
    cpdef Error(self):
        return dvectond(self.sol[0].sys[0].Error())
    cpdef GetRes(self):
        return dvectond(self.sol[0].sys[0].GetRes())
    cpdef DoIter(self,d relaxfac):
        self.sol[0].DoIter(relaxfac)
    cpdef SetRes(self,n.ndarray[n.float64_t,ndim=1] res):
        self.sol[0].sys[0].SetRes(dndtovec(res))
    cpdef GetResVar(self,_type,freqnr):
        if _type=='pres':
            return dvectond(self.tube[0].GetResAt(3,freqnr))
        elif _type=='rho':
            return dvectond(self.tube[0].GetResAt(0,freqnr))
        elif _type=='temp':
            return dvectond(self.tube[0].GetResAt(2,freqnr))
        elif _type=='volu':
            return dvectond(self.tube[0].GetResAt(1,freqnr))
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

