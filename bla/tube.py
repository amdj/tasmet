#!/usr/bin/python

from cmath import *
import numpy as n
from copy import *
from gridpoint import *

class tube:
	def __init__(s,T,Nf=1):
		s.gp=[]
		s.Nf=Nf
		s.Ns=2*Nf+1
		s.T=T
		for i in range(10):
			s.gp.append(gridpoint(s.Nf,1.))
	def __repr__(s):
		return "Tube class"
# 	def d_rho_u_dx(s,i):
# 		rho_u_ip1=
	
	
