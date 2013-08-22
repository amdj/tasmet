#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 21:05:46 2012

@author: anne
"""
import numpy as n
from math import *
class air:
	def __init__(s,T0=293.15):
		s.Rs=287.
		s.T0=T0
		s.cp0c=1047.63657
		s.cp1c=-0.372589265
		s.cp2c=9.45304214E-4
		s.cp3c=-6.02409443E-7
		s.cp4c=1.2858961E-10
		s.cp0=s.cp_T(T0)
		s.h0=s.h_T(T0)
	def kappa_T(self,T):
# 		print "Temporarily set kappa to constant in materials.py"
		return -0.00227583562+1.15480022E-4*T**1-7.90252856E-8*T**2+4.11702505E-11*T**3-7.43864331E-15*T**4
# 			return 0.0257*T**0
	def rho(self,p,T):
# 		print "Materials15,T=",T
# 			raise IOError('Zero temperature')
# 		print "materials15,T=",T
		return p/self.Rs/T
	def p_rho_T(self,rho,T):
		return rho*self.Rs*T
	def mu_T(self,T):
		return -8.38278E-7+8.35717342E-8*T**1-7.69429583E-11*T**2+4.6437266E-14*T**3-1.06585607E-17*T**4
	def cp_T(self,T):
		return self.cp0c+self.cp1c*T**1+self.cp2c*T**2+self.cp3c*T**3+self.cp4c*T**4
	def intcp_over_TdT(self,T0,T):
		return self.cp0c*n.log(T/T0)+self.cp1c*(T-T0)+0.5*self.cp2c*(T**2-T0**2)+self.cp3c*(T**3-T0**3)/3.+self.cp4c*(T**4-T0**4)/4.
	def cv_T(s,T):
		return s.cp_T(T)-s.Rs
	def h_T(s,T):
		#We have to integrate the cp_T function, for this, it is made a bit easier
		a=(s.cp0c*(T)+0.5*s.cp1c*(T**2)+
				s.cp2c*(1/3.0)*(T)**3+
				s.cp3c*0.25*(T)**4+
				s.cp4c*(0.2)*(T)**5)
# 		a=1004.5*T
		return a
	def e_T(s,T):
		return s.h_T(T)-s.Rs*(T)
# 		return (1004.5-287)*T
	def gamma_T(self,T):
		cp=self.cp_T(T)
		return cp/(cp-self.Rs)
	def cm_T(self,T):
# 		 print "materials24,T equals:",T
# 		if type(T)!=type(1.0):
# 			print T
# 			for t in range(T):
# 				if T[t]<0:
# 					raise IOError('Temp < 0!')

		cm=n.sqrt(self.gamma_T(T)*self.Rs*T)
#		 print "materials26,cm equals:",cm
		return cm
