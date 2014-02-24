'''
Created on Aug 5, 2013

@author: anne
'''
import numpy as n
class testbc_constantcs:
	def testInit(self):
		self.pm=101325.
		self.Tm=293.15
		self.testp=n.zeros((1,self.Nf+1),complex)
		self.testT=n.zeros((1,self.Nf+1),complex)
		self.testT[0,0]=self.Tm
		self.testp[0,0]=self.pm
		self.testU=n.zeros((1,self.Nf+1),complex)
		self.upiston=n.zeros((1,self.Nf+1),complex)
		self.upiston[0,1]=1.
		self.rhom=self.pm/287/self.Tm
		self.Vpiston=n.zeros(self.upiston.shape,complex)
		self.Vpiston[0,1:]=self.upiston[0,1:]/1j/self.omega[0,1:]
		self.Vpiston[0,0]=self.x[1]
	def leftBc3(self,error):
		i=0
		rhoblk,Ublk,pblk,Tblk=self.blocks(i)
		Vpiston_td=self.idft1long(self.Vpiston)
		rhoudotn=(self.mf()[i]+self.mf()[i+1])
		dmdt=1j*self.omega*self.dft(self.rho.getTdata()[i]*Vpiston_td)
		#Density equation
		
		error[rhoblk[0]:rhoblk[1]]=dmdt+rhoudotn
		error[rhoblk[0]]=self.rho.getAdata()[0,0]-self.rhom

		#Momentum balance
		dmfdt=1j*self.omega*self.dft(self.U.getTdata()[i]*Vpiston_td)
		mominleft_td=self.rho.getTdata()[i]*Vpiston_td**2
		mominleft=self.dft(mominleft_td)
		pdotn=-self.p()[i]+0.5*(self.p()[i]+self.p()[i+1])
		rhouudotn=0.5*(self.momf()[i]+self.momf()[i+1])-mominleft

		error[Ublk[0]:Ublk[1]]=dmfdt+pdotn#+rhouudotn
 		error[Ublk[0]]=self.mf()[i,0]

		#dEtotdt=1j*self.omega*self.Etot()[i]*Vi
		#Hleft_td=minleft_td*(self.m.h_T(self.T.getTdata()[i])+0.5*uleft_td*uleft_td)
		#Hleft=self.dft(Hleft_td)

		#Couple pressure, density and temperature by adiabatic state equations
		#Hleft_td=
		#Terr_td=1/self.m.Rs*self.m.intcp_over_TdT(self.Tm, self.T.getTdata()[i])-self.p.getTdata()[i]/101325.
		#error[Tblk[0]:Tblk[1]]=self.dft(Terr_td)
		#EASIER: isentropic law!
		error[Tblk[0]:Tblk[1]]=self.T()[i]-self.Tm*self.dft((self.rho.getTdata()[i]/self.rhom)**0.4)
		error[Tblk[0]]=self.T()[i,0]-293.15
		#Perfect gas law equation of state
		error[pblk[0]:pblk[1]]=self.p()[i]-self.dft(self.m.p_rho_T(self.rho.getTdata()[i], self.T.getTdata()[i]))

		return error
	def leftBclin(self,error):
		i=0
		rhoblk,Ublk,pblk,Tblk=self.blocks(i)
		Vpiston_td=self.idft1long(self.Vpiston)
		drhoudxlin=self.rhom*(self.U()[i+1]-self.U()[i])/(self.x[i+1]-self.x[i])
		dmdt=1j*self.omega*self.rho()[i]
		#Density equation
		
		error[rhoblk[0]:rhoblk[1]]=dmdt+drhoudxlin
		error[rhoblk[0]]=self.rho.getAdata()[0,0]-self.rhom

		#Momentum balance

		error[Ublk[0]:Ublk[1]]=self.U()[i]-self.upiston
		#EASIER: isentropic law!
		error[Tblk[0]:Tblk[1]]=self.T()[i]-self.Tm*self.dft((self.rho.getTdata()[i]/self.rhom)**0.4)
		error[Tblk[0]]=self.T()[i,0]-self.Tm
		#Perfect gas law equation of state
		error[pblk[0]:pblk[1]]=self.p()[i]-self.dft(self.m.p_rho_T(self.rho.getTdata()[i], self.T.getTdata()[i]))

		return error
		
	def rightBc(self,error):
		i=self.gridpoints-1
		rhoblk,Ublk,pblk,Tblk=self.blocks(i)
		#Implementation of continuity equation
		Vi=self.x[i]-self.x[i-1]
		dmdt=1j*self.omega*self.rho()[i]*Vi
		rhoudotn=-0.5*(self.mf()[i-1]+self.mf()[i])
		error[rhoblk[0]:rhoblk[1]]=dmdt+rhoudotn

		#Implementation of momentum equation
		dmfdt=1j*self.omega*self.mf()[i]*Vi
		rhouudotn=-0.5*(self.momf()[i-1]+self.momf()[i])
		pdotn=self.p()[i]-0.5*(self.p()[i]+self.p()[i-1])
		#FIXTHIS
		error[Ublk[0]:Ublk[1]]=dmfdt+pdotn#+rhouudotn
		#error[Ublk[0]]=self.momf()[i,0]

		#Implementation of energy equation
		#dEtotdt=1j*self.omega*self.Etot()[i]*Vi
		#error[Tblk[0]:Tblk[1]]=dEtotdt-0.5*(self.Hf()[i]+self.Hf()[i-1])
		#Implementation isentropic law
		error[Tblk[0]:Tblk[1]]=self.T()[i]-self.Tm*self.dft((self.rho.getTdata()[i]/self.rhom)**0.4)
		error[Tblk[0]]=self.T()[i,0]-self.Tm
		#Implementation of equation of state
		error[pblk[0]:pblk[1]]=self.p()[i]-self.dft(self.m.p_rho_T(self.rho.getTdata()[i], self.T.getTdata()[i]))
		return error
class testing:
	pass
