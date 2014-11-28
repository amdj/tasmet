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
    def getFreq(self):
        return self.tube.getFreq()
    def getNf(self):
        return self.tube.getNf()
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
    def __init__(self,segs):
        self.segs=segs
    def getx(self):
        try:
            return self.x
        except:
            L=0
            x=[]
            for seg in self.segs:
                x=n.concatenate((x,seg.getx()+L))
                L+=seg.L()
            self.x=x
            return x
    def L(self):
        return self.x[-1]
    def getNf(self):
        return self.segs[0].getNf()
    def getFreq(self):
        return self.segs[0].getFreq()
    def quantityAtTimeAndPlaceInterp(self,t,x,quant="pres"):
        self.getx()
        assert(x>0 and x<self.L())
        iR=0
        while self.x[iR]<x:
            iR+=1
        iL=iR-1
        xL=self.x[iL]
        xR=self.x[iR]        
        relL=(x-xL)/(xR-xL)
        return self.quantityAtTimeAndPlace(t,iL,quant)*relL+(1-relL)*self.quantityAtTimeAndPlace(t,iR,quant)            
        
    def quantityAtTimeAndPlace(self,t,i,quant="pres"):

        timeval=0.0
        omg=2*n.pi*self.getFreq()
        for j in range(self.getNf()+1):
            timeval+=(fouriercoefs[j]*n.exp(j*1j*omg*t)).real
        return timeval

    def fouriercoefs(self,quant):
        if quant=="pres":
            try:
                return self.prescoefs
            else:
                self.prescoefs=fouriercoefs(quant)
                return self.prescoefs
        elif quant=="velo":
            try:
                return self.velocoefs
            else:
                self.velocoefs=fouriercoefs(quant)
                return self.velocoefs
        elif quant=="volu":
            try:
                return self.volucoefs
            else:
                self.volucoefs=fouriercoefs(quant)
                return self.volucoefs
        
            
    def getfouriercoefs(self,quant):
        fouriercoefs=[]
        for freqnr in range(self.getNf()+1):
            fouriercoefs.append(self.getvar(freqnr,quant)[i])
        return fouriercoefs

    def quantityAtTime(self,t,quant="pres"):
        fouriercoefs=[]
        result=n.zeros(self.getx().size)
        omg=2*n.pi*self.getFreq()
        for freqnr in range(self.getNf()+1):
            fouriercoefs.append(self.getvar(freqnr,quant))
        for freqnr in range(self.getNf()+1):
            result+=(fouriercoefs[freqnr]*n.exp(1j*freqnr*omg*t)).real
        return result
    def dragCoef(self,freqnr):
        dc=[]    
        for seg in self.segs:
            dc=n.concatenate((dc,seg.dragCoef(freqnr)))
        return dc
    def getvar(self,freqnr,name="pres"):
        if name=="volu":
            return self.getU(freqnr)
        var=[]
        for seg in self.segs:
            var=n.concatenate((var,seg.getvar(freqnr,name)))
        return var
    def getSf(self):
        Sf=[]    
        for seg in self.segs:
            Sf=n.concatenate((Sf,seg.getSf()))
        return Sf
    def getS(self):
        S=[]    
        for seg in self.segs:
            S=n.concatenate((S,seg.getS()))
        return S
    def getrh(self):
        rh=[]    
        for seg in self.segs:
            rh=n.concatenate((rh,seg.getrh()))
        return rh
    def getphi(self):
        phi=[]    
        for seg in self.segs:
            phi=n.concatenate((phi,seg.getphi()))
        return phi
        
    def getp(self,i):
        return self.getvar(i,"pres")
    def Edot(self,i):
        return (0.5*self.getp(i)*self.getU(i).conjugate()).real
    def getrho(self,i):
        rho=[]
        for seg in self.segs:
            rho=n.concatenate((rho,seg.getrho(i)))
        return rho
    def getT(self,i):
        T=[]
        for seg in self.segs:
            T=n.concatenate((T,seg.getT(i)))
        return T
    def getTs(self,i):
        Ts=[]
        for seg in self.segs:
            Ts=n.concatenate((Ts,seg.getTs(i)))
        return Ts
    def getHtot(self):
        Htot=[]
        for seg in self.segs:
            Htot=n.concatenate((Htot,seg.getHtot()))
        return Htot

    def getU(self,i):
        try:
            return self.U[i]
        except:
            print("Making first U")
            self.U=[]
            for i in range(self.getNf()+1):
                self.U.append([])
                for seg in self.segs:
                    self.U[i]=n.concatenate((self.U[i],seg.getU(i)))
        return self.U[i]
    def getu(self,i):
