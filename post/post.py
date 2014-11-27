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
        self.tube=tube
        self.nr=nr
    def getx(self):
        return self.tube.getx(self.nr)
    def L(self):
        return self.getx()[-1]
    def getSf(self):
        return self.tube.getSf(self.nr)
    def getS(self):
        return self.tube.getS(self.nr)
    def getphi(self):
        return self.tube.getphi(self.nr)
    def getrh(self):
        return self.tube.getrh(self.nr)
    def dragCoef(self,freqnr):
        if(freqnr==0):
            return self.tube.dragCoef(0,self.nr)
        else:
            return self.tube.dragCoef(2*freqnr-1,self.nr)+1j*self.tube.dragCoef(2*freqnr,self.nr)
    def getvar(self,i,name):
        if(i==0):
            return self.tube.getResVar(name,i,self.nr)
        else: 
            return self.tube.getResVar(name,2*i-1,self.nr)+1j*self.tube.getResVar(name,2*i,self.nr)       
    def getp(self,i):
        return self.getvar(i,"pres")
    def getrho(self,i):
        return self.getvar(i,"rho")
    def getT(self,i):
        return self.getvar(i,"temp")
    def getTs(self,i):
        return self.getvar(i,"stemp")
    def getU(self,i):
        return self.getvar(i,"volu")
    def getHtot(self):
        return self.tube.getHtot(self.nr)
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
        self.systems=systems
    def getx(self):
        x=[]
        L=0
        for sys in self.systems:
            x=n.concatenate((x,sys.getx()+L))
            L+=sys.L()
        self.x=x
        return x
        
    def dragCoef(self,freqnr):
        dc=[]    
        for sys in self.systems:
            dc=n.concatenate((dc,sys.dragCoef(freqnr)))
        return dc
    def getSf(self):
        Sf=[]    
        for sys in self.systems:
            Sf=n.concatenate((Sf,sys.getSf()))
        return Sf
    def getS(self):
        S=[]    
        for sys in self.systems:
            S=n.concatenate((S,sys.getS()))
        return S
    def getrh(self):
        rh=[]    
        for sys in self.systems:
            rh=n.concatenate((rh,sys.getrh()))
        return rh
    def getphi(self):
        phi=[]    
        for sys in self.systems:
            phi=n.concatenate((phi,sys.getphi()))
        return phi
        
    def getp(self,i):
        p=[]
        for sys in self.systems:
            p=n.concatenate((p,sys.getp(i)))
        return p
    def Edot(self,i):
        return (0.5*self.getp(i)*self.getU(i).conjugate()).real
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
