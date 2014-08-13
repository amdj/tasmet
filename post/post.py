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
    def __init__(self,tube):
        nlpost.__init__(self,tube.getFreq(),tube.getNf())

        self.rhoi=[]
        self.p=[]
        self.T=[]
        self.Ts=[]
        self.rho=[]
        self.U=[]
        self.x=tube.getx()

        self.Ts.append(tube.getResVar('stemp',0))
        self.T.append(tube.getResVar('temp',0))
        self.p.append(tube.getResVar('pres',0))
        self.rho.append(tube.getResVar('rho',0))
        self.U.append(tube.getResVar('volu',0))
        #Compute Fubini solution
        Nf=self.Nf
        for i in range(1,Nf+1):
            self.rho.append(tube.getResVar("rho",2*i-1)+1j*tube.getResVar("rho",2*i))
            self.p.append(tube.getResVar("pres",2*i-1)+1j*tube.getResVar("pres",2*i))
            self.T.append(tube.getResVar("temp",2*i-1)+1j*tube.getResVar("temp",2*i))
            self.U.append(tube.getResVar("volu",2*i-1)+1j*tube.getResVar("volu",2*i))
            self.Ts.append(tube.getResVar("stemp",2*i-1)+1j*tube.getResVar("stemp",2*i))
        # self.iDFT[:,0]=1.
    def getp(self,i):
        return self.p[i]
    def getrho(self,i):
        return self.rho[i]
    def getT(self,i):
        return self.T[i]
    def getTs(self,i):
        return self.Ts[i]
    def getU(self,i):
        return self.U[i]

