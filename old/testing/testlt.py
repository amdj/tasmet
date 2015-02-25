#!/usr/bin/python

from tmtubes.linear import * 
from lintube_python import *
import numpy as n
Nf=1
gp=10
Up=1.
freq=60
t=lintube_py(gp,Nf,freq,Up)
solinit=t.getResult()
initEr=t.getError()
# print solinit
f=100.
c0=sqrt(1.4*287*293.15)
lambda_=c0/f
t.setResult(solinit)
sol=new_raph_numpy(t.Solve,solinit,1e-6,verbose=True)
#p=t.getpResult(1)+1j*t.getpResult(2)
#U=t.getUResult(1)+1j*t.getUResult(2)
#x=t.getx()
