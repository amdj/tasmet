#!/usr/bin/python

# This file provides the class TaSMETModel
import pickle
from os.path import isfile
import TASMET as tasmet

class TaSMETModel(object):
    def __init__(self):
        self.Nf=6
        self.freq=100.
        self.T0=293.15
        self.p0=101325.
        self.gas='air'
        self.gases={'air','helium','nitrogen'}
        self.m0=0.
        self.EngineSystem=False

        self.segs=[]
        self.cons=[]

    @classmethod
    def load(self,filename):
        if(isfile(filename)):
            return pickle.load(open(filename,'rb'))
        else:
            print('%s not found.')
            
    def save(self,filename):
        f=open(filename,'wb')
        pickle.dump(self,f)

    def buildModel(self):
        gc=tasmet.Globalconf(self.Nf,self.freq,self.gas,self.T0,self.p0)
        if self.EngineSystem:
            sys=tasmet.EngineSystem(gc)
        else:
            sys=tasmet.TaSystem(gc)

        for seg in self.segs:
            sys.addSeg(self.makeSeg(seg))
        return sys

if __name__=='__main__':
    t=TasMETModel()
    t.save('test.dat')
