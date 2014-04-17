#!/usr/bin/python

import numpy as n

class nonlpost:
    def __init__(self,freq,Nf):
        self.freq=freq
        self.Nf=Nf
        self.Ns=2*Nf+1
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
        # self.iDFT[:,0]=1.


