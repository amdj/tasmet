#!/usr/bin/python

import numpy as n

class nlpost(object):
    def __init__(self,freq,Nf):
        self.Nf=Nf
        self.freq=freq
        self.Ns=2*self.Nf+1
        self.mkfDFT()
        self.mkiDFT()

    def mkfDFT(self):
        Ns=self.Ns
        self.fDFT=n.zeros((Ns,Ns),float)
        self.fDFT[0,:]=1.0/Ns
        for i in range(1,self.Nf+1):
            for j in range(self.Ns):
                self.fDFT[2*i-1,j]=2.0*n.cos(2.0*n.pi*i*j/self.Ns)/Ns
                self.fDFT[2*i,j]=-2.0*n.sin(2.0*n.pi*i*j/self.Ns)/Ns
    def mkiDFT(self):
        Ns=self.Ns
        Nf=self.Nf
        self.iDFT=n.zeros((Ns,Ns),float)
        self.iDFT[:,0]=1.0
        for k in range(0,Ns):
            for r in range(1,Nf+1):
                self.iDFT[k,2*r-1]=n.cos(2.0*n.pi*r*k/Ns)
                self.iDFT[k,2*r]=-n.sin(2.0*n.pi*r*k/Ns)                
    

class nonlinpost(nlpost):
    def __init__(self,tube,nr=0):
        nlpost.__init__(self,tube.getFreq(),tube.getNf())
        self.ntubes=tube.nTubes()
        self.rhoi=[]
        self.p=[]
        self.T=[]
        self.Ts=[]
        self.rho=[]
        self.U=[]
        self.x=tube.getx(nr)
        self.Sf=tube.getSf(nr)
        self.L=tube.getL(nr)
        self.Ts.append(tube.getResVar('stemp',0,nr))
        self.T.append(tube.getResVar('temp',0,nr))
        self.p.append(tube.getResVar('pres',0,nr))
        self.rho.append(tube.getResVar('rho',0,nr))
        self.U.append(tube.getResVar('volu',0,nr))
        #Compute Fubini solution
        self.Htot=tube.getHtot(nr)
        Nf=self.Nf
        for i in range(1,Nf+1):
            self.rho.append(tube.getResVar("rho",2*i-1,nr)+1j*tube.getResVar("rho",2*i,nr))
            self.p.append(tube.getResVar("pres",2*i-1,nr)+1j*tube.getResVar("pres",2*i,nr))
            self.T.append(tube.getResVar("temp",2*i-1,nr)+1j*tube.getResVar("temp",2*i,nr))
            self.U.append(tube.getResVar("volu",2*i-1,nr)+1j*tube.getResVar("volu",2*i,nr))
            self.Ts.append(tube.getResVar("stemp",2*i-1,nr)+1j*tube.getResVar("stemp",2*i,nr))
        # self.iDFT[:,0]=1.
    def getp(self,i):
        return self.p[i]
    def getSf(self):
        return self.Sf
    def getrho(self,i):
        return self.rho[i]
    def getT(self,i):
        return self.T[i]
    def getTs(self,i):
        return self.Ts[i]
    def getU(self,i):
        return self.U[i]
    def getHtot(self):
        return self.Htot            
    def pressureprofile(self,pos,Nperiods=2,ns=100):
        p=[]
        for k in range(self.Nf+1):
            p.append(self.getp(k)[pos])
        t=n.linspace(0,Nperiods/self.freq,ns)
        pt=p[0]*n.ones(t.shape)
        omg=2*n.pi*self.freq
        for k in range(1,self.Nf+1):
            pt+=(p[k]*n.exp(1j*k*(omg*t))).real        
        return (t,pt)

class combinedsys:
    def __init__(self,systems):
        x=[]
        L=0
        self.systems=systems
        for sys in systems:
            x=n.concatenate((x,sys.x+L))
            L+=sys.L

        self.x=x
    def getSf(self):
        Sf=[]    
        for sys in self.systems:
            Sf=n.concatenate((Sf,sys.getSf()))
        return Sf

        
    def getp(self,i):
        p=[]
        for sys in self.systems:
            p=n.concatenate((p,sys.getp(i)))
        return p
    def getp(self,i):
        p=[]
        for sys in self.systems:
            p=n.concatenate((p,sys.getp(i)))
        return p
    def getrho(self,i):
        rho=[]
        for sys in self.systems:
            rho=n.concatenate((rho,sys.getrho(i)))
        return rho
    def getT(self,i):
        T=[]
        for sys in self.systems:
            T=n.concatenate((T,sys.getT(i)))
        return T
    def getTs(self,i):
        Ts=[]
        for sys in self.systems:
            Ts=n.concatenate((Ts,sys.getTs(i)))
        return Ts
    def getHtot(self):
        Htot=[]
        for sys in self.systems:
            Htot=n.concatenate((Htot,sys.getHtot()))
        return Htot

    def getU(self,i):
        U=[]
        for sys in self.systems:
            U=n.concatenate((U,sys.getU(i)))
        return U
