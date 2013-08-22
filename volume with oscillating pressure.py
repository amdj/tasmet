# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>
from numpy import *
from numpy.linalg import *
from matplotlib.pyplot import *
#---------------------------------------------------------- from scipy. import *
from volume import *
from basefuncs import *
Nf=15
v=volume(1.,Nf=Nf)
V=0.01
S=0.1
v.Init(V=V,S=S)

omg=2*pi
R=287.

pm=101325.

Tm=293.15
rhom=pm/Tm/R
p1=80000.
gam=1.4
v.addbc('pres',[pm,p1])

rho1=p1/(gam*R*Tm)
T1=Tm*(p1/pm-rho1/rhom)
u1=-V*1j*omg*rho1/(rhom*S)
#guess1=n.array([pm,p1,0,0,rhom,rho1,0,0,Tm,T1,0,0.,0,u1,0,0])
# 7 harmonics
guess=n.array([pm,p1,0,0,0,0,0,0,
						1e-6,u1,0,0,0,0,0,0,
							Tm,T1,0,0,0,0,0,0,
					rhom,rho1,0,0,0,0,0,0])
# 9 harmonics
guess=n.array([pm,p1,0,0,0,0,0,0,0,0,
						1e-6,u1,0,0,0,0,0,0,0,0,
							Tm,T1,0,0,0,0,0,0,0,0,
					rhom,rho1,0,0,0,0,0,0,0,0])
# 11 harmonics
guess=n.array([pm,p1,0,0,0,0,0,0,0,0,0,0,
						1e-6,u1,0,0,0,0,0,0,0,0,0,0,
							Tm,T1,0,0,0,0,0,0,0,0,0,0,
					rhom,rho1,0,0,0,0,0,0,0,0,0,0])
# 15 harmonics
guess=n.array([pm,p1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
						1e-6,u1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
							Tm,T1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
					rhom,rho1,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

guess=n.array([pm,p1,0,0,0,0,1e-6,u1,0,0,0,0,Tm,T1,0,0,0,0,rhom,rho1,0,0,0,0])
del guess
guess=n.array([pm,0 ,0,0,0,0,0   ,0 ,0,0,0,0,0 ,0 ,0,0,0,0,rhom,0   ,0,0,0,0])
#guess=n.array([pm,p1,0.0001,u1,Tm,T1,rhom,rho1])
#guess=1e-6*n.ones(((Nf+1)*4,),complex)
E=v.err(guess)

#kl=1e-6
#ku=-5*kl
#kt=-1*kl
#kr=-5*kl
def complextofull(complexval):
		length=complexval.shape[0]

		full=n.zeros((2*length,),float)
		full[0:2*length:2]=complexval.real[:]
		full[1:2*length+1:2]=complexval.imag[:]
		return full

def fulltocomplex(full):
		length=full.shape[0]/2
		complexval=n.zeros((length,),complex)

		complexval[0:length]=full[0:2*length:2]
		complexval[0:length]+=1j*full[1:2*length+1:2]
		return complexval

#print E[6:8]

#du1=ku*E[Nf+1:2*(Nf+1)]/rhom/S
#dT1=kt*E[2*(Nf+1):3*(Nf+1)]/rhom/S/1000
#drho1=kr*E[3*(Nf+1):4*(Nf+1)]
#v.phys.u.setAdata(v.phys.u()+du1)
#v.phys.temp.setAdata(v.phys.temp()+dT1)
#v.phys.rho.setAdata(v.phys.rho()+drho1)

print norm(E)

# <codecell>

sol=complextofull(guess)

def error(guessreal):
		complexguess=fulltocomplex(guessreal)
		error=v.err(complexguess)
		errorreal=complextofull(error)
		return errorreal
from scipy.optimize import *
#sol=new_raph_numpy(error,sol,1e-6,1e-6,maxiter=100)
#sol=newton_krylov(error,sol,verbose=False)
# sol=anderson(error,sol,verbose=False)
guess=fulltocomplex(sol)

# <codecell>

time=v.phys.time.getProfile(400,4)
pres=n.asarray(v.phys.pres.getProfile(400,4))
temp=n.asarray(v.phys.temp.getProfile(400,4))
figure(num=None, figsize=(15, 6), dpi=80, facecolor='w', edgecolor='k')
plot(time,temp)

tempexact=Tm*(pres/pm)**((gam-1)/(gam))
tempsin=293.15+(T1*n.exp(1j*omg*time)).real
hold('on')

plot(time,tempexact,'red')
plot(time,tempsin,'green')

legend(('Numerical solution','Exact solution','Linearized exact'),'lower left')

# <codecell>

print tempexact[25]
print sum(tempexact)/tempexact.shape[0]
print Tm*(pres[25]/pm)**((gam-1)/gam)

# <codecell>

u=v.phys.u.getProfile(100,2)
t=v.phys.time.getProfile(100,2)
plot(t,u)

# <codecell>

plot((abs(v.phys.u()[:])))

# <codecell>

plot(v.phys.pres.getProfile(100,3))

# <codecell>


