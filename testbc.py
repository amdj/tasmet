'''
Created on Aug 5, 2013

@author: anne
'''
import numpy as n
class testbc:
	def leftBc3(self,error):
		i=0
		rhoblk,Ublk,pblk,Tblk=self.blocks(i)
		Vpiston=self.upiston*self.Sf[0]
		
		Vpiston_td=self.idft1long(self.upiston*self.Sf[0])
		rhoudotn=0.5*(self.mf()[i]+self.mf()[i+1])
		dmdt=1j*self.omega*self.dft(self.rho.getTdata()[i]*Vpiston_td)
		#Density equation
		
		error[rhoblk[0]:rhoblk[1]]=dmdt+rhoudotn
		error[rhoblk[0]]=self.rho.getAdata()[0,0]-1.2

		#Momentum balance
		dmfdt=1j*self.omega*self.dft(self.mf.getTdata()[i]*Vpiston_td)
		mominleft_td=self.rho.getTdata()[i]*Vpiston_td**2/self.SfT[i]
		mominleft=self.dft(mominleft_td)
		pdotn=-self.Sf[i]*self.p()[i]+0.5*(self.Sf[i]*self.p()[i]+self.Sf[i+1]*self.p()[i+1])
		rhouudotn=0.5*(self.momf()[i]+self.momf()[i+1])-mominleft+pdotn
		#No p dot n
		#HIER ZIT HET NOG NIET HELEMAAL GOED??
		#FIXTHIS
		error[Ublk[0]:Ublk[1]]=dmfdt+rhouudotn+pdotn
 		error[Ublk[0]]=self.mf()[i,0]

		#dEtotdt=1j*self.omega*self.Etot()[i]*Vi
		#Hleft_td=minleft_td*(self.m.h_T(self.T.getTdata()[i])+0.5*uleft_td*uleft_td)
		#Hleft=self.dft(Hleft_td)

		#Couple pressure, density and temperature by adiabatic state equations
		#Hleft_td=
		Terr_td=1/self.m.Rs*self.m.intcp_over_TdT(293.15, self.T.getTdata()[i])-self.p.getTdata()[i]/101325.
		error[Tblk[0]:Tblk[1]]=self.dft(Terr_td)
		error[Tblk[0]]=self.T()[i,0]-293.15
		#Perfect gas law equation of state
		error[pblk[0]:pblk[1]]=self.p()[i]-self.dft(self.m.p_rho_T(self.rho.getTdata()[i], self.T.getTdata()[i]))

		return error

	def rightBc(self,error):
		i=self.gridpoints-1
		rhoblk,Ublk,pblk,Tblk=self.blocks(i)
		#Implementation of continuity equation
		Vi=self.Sf[i]*(self.x[i]-self.x[i-1])
		dmdt=1j*self.omega*self.rho()[i]*Vi
		rhoudotn=-0.5*(self.mf()[i-1]+self.mf()[i])
		error[rhoblk[0]:rhoblk[1]]=dmdt+rhoudotn

		#Implementation of momentum equation
		dmfdt=1j*self.omega*self.mf()[i]*Vi
		rhouudotn=-0.5*(self.momf()[i-1]+self.momf()[i])
		pdotn=self.Sf[i]*self.p()[i]-0.5*(self.Sf[i]*self.p()[i]+self.Sf[i-1]*self.p()[i-1])
		#FIXTHIS
		error[Ublk[0]:Ublk[1]]=dmfdt+pdotn+rhouudotn
		#error[Ublk[0]]=self.momf()[i,0]

		#Implementation of energy equation
		dEtotdt=1j*self.omega*self.Etot()[i]*Vi
		error[Tblk[0]:Tblk[1]]=dEtotdt-0.5*(self.Hf()[i]+self.Hf()[i-1])
		error[Tblk[0]]=self.T()[i,0]-293.15
		#Implementation of equation of state
		error[pblk[0]:pblk[1]]=self.p()[i]-self.dft(self.m.p_rho_T(self.rho.getTdata()[i], self.T.getTdata()[i]))

		return error
class testing:
	def testInit(self):
		self.testp=n.zeros((1,self.Nf+1),complex)
		self.testT=n.zeros((1,self.Nf+1),complex)
		self.testT[0,0]=293.15
		self.testp[0,0]=101325
		#self.testp[0,1]=0
		self.testU=n.zeros((1,self.Nf+1),complex)
		self.upiston=n.zeros((1,self.Nf+1),complex)
		self.upiston[0,1]=1e-3*1j/200
		#self.testT[0]=293.15
