#!/usr/bin/python

from cmath import *
import mpmath as mp
import math as m
import numpy as n

def Zp(omega,S,rho0,c0):
	D=2*m.sqrt(S/pi)
	k=omega/c0
	kD=k*D
# 	print "basefuncs.py8, k*D=",kD
	R1=1-2*mp.besselj(1,kD)/(kD)
	X1=2*mp.struveh(1,kD)/(kD)
	Z0=rho0*c0
	Zp=Z0*(R1+1j*X1)
	return Zp

# omega=500*2*pi+2j
# S=pi*(.5*2.1e-2)**2
# rho0=1.2
# c0=343

# print Zp(omega,S,rho0,c0)

def broyden(fun,guess,guess2=None,funtol=1e-3,reltol=1e-1,maxiter=10,verbose=False):
	# Expected, input a float numpy matrix guess vector of shape (L,1)
	# Expected: tol as a float
	if verbose:
		print "Broyden iteration on function ", fun ,"started"
	nloop=0
	nconverg=0
	sol=guess
	nsys=sol.shape[0]
	fx=fun(n.copy(sol))
	x1=sol.copy()
	J=n.zeros((nsys,nsys),float)
	normx=n.linalg.norm(sol)
	if guess2 is None:
		for j in xrange(nsys):
			x2=n.copy(x1)
			unitj=n.zeros((nsys,))
			unitj[j]=0.01*reltol*normx
			x2=x2+unitj
			fx2=fun(x2)
			dfdxj=(fx2-fx)/(x2[j]-x1[j])
			J[:,j]=dfdxj.ravel()
		if n.linalg.det(J)==0:
#			print n.diag(J)
			print "Jabobian equals zero for function" +  str(fun) + ". Exiting solving procedure."
			return sol,0,0
			raise IOError('Determinant of Jacobian becomes zero')
		dx=(-n.linalg.solve(J,fx)).ravel()
		sol=sol+dx
	else:
		pass
# 		fx2=fun(guess2)
# 		for i in xrange(nsys):
# 			diagJ=(fx2-fx)
# 			J[i,i]
	relerror=n.linalg.norm(dx)/n.linalg.norm(sol)
	funerror=n.linalg.norm(fx)
	if relerror < reltol and funerror<funtol:
		nloop=maxiter
###########
	while nloop<maxiter:
		fxold=fx.copy()
		fx=fun(sol)
		df=fx-fxold
		denom=n.linalg.norm(dx)**2
		J=J+(df-J*dx)*dx.T/denom
		dx=-n.linalg.solve(J,fx)
		sol=sol+dx
		error=n.asarray(dx)
		nloop+=1
		nconverg+=1
		if nloop==maxiter-1:
			print "Warning: maximum number of iterations reached, but convergence is not reached!"
		relerror=n.linalg.norm(dx)/n.linalg.norm(sol)
		funerror=n.linalg.norm(fx)
		if verbose==True:
			print "Iteration: %g. Relative error: %.2e. Function error: %.2e" %(nloop,relerror,funerror)
# 		print "Iteration ,",nloop,", current relative error:",relerror, ", current func. value norm: ",funnorm
		if relerror<reltol and funerror<funtol:
			nloop=maxiter

	if verbose:
		print "Iteration on function ", fun ,"converged"
	return sol,nconverg,relerror
############################################################

def broyden_twoguess(fun,guess1,guess2,funtol=1e-3,reltol=1e-1,maxiter=10,verbose=False):
	# Expected, input a float numpy matrix guess vector of shape (L,1)
	# Expected: tol as a float
	nloop=0
	nconverg=0
	sol=n.mat(guess1)
	nsys=len(sol)
	fx=fun(n.array(sol))
	sol2=n.mat(guess2)
	nsys=len(sol)
	fx2=fun(n.array(sol2))
	x1=guess1
	J=n.mat(n.zeros((nsys,nsys),float))

	normx=n.linalg.norm(sol)
	for j in xrange(nsys):
		x2=x1.copy()
		unitj=n.zeros((nsys,1))
		unitj[j]=0.01*reltol*normx
		unitj=n.mat(unitj)
#  			print unitj
		x2=x2+unitj
		fx2=fun(n.array(x2))
		dfdxj=(fx2-fx)/(x2[j]-x1[j])
		J[:,j]=dfdxj
	if n.linalg.det(J)==0:
		raise IOError('Determinant of Jacobian becomes zero')
	dx=-n.linalg.solve(J,fx)
	sol=sol+dx
	relerror=n.linalg.norm(dx)/n.linalg.norm(sol)
	funerror=n.linalg.norm(fx)
	if relerror < reltol and funerror<funtol:
		nloop=maxiter
###########
	if verbose:
		print('Broyden iteration on function started')
	while nloop<maxiter:
		fxold=fx.copy()
		fx=fun(sol)
		df=fx-fxold
		denom=n.linalg.norm(dx)**2
		J=J+(df-J*dx)*dx.T/denom
		dx=-n.linalg.solve(J,fx)
		sol=sol+dx
		error=n.asarray(dx)
		nloop+=1
		nconverg+=1
		if nloop==maxiter-1:
			print 'Warning: maximum number of iterations reached, but convergence is not reached!'
		relerror=n.linalg.norm(dx)/n.linalg.norm(sol)
		funerror=n.linalg.norm(fx)
		if verbose==True:
			print "Iteration: %g. Relative error: %.2e. Function error: %.2e" %(nloop,relerror,funerror)
# 		print "Iteration ,",nloop,", current relative error:",relerror, ", current func. value norm: ",funnorm
		if relerror<reltol and funerror<funtol:
			nloop=maxiter

	if verbose:
		print "Iteration on function ", fun ,"converged"
	return n.array(sol),nconverg,relerror
############################################################
####################################################################
def broyden_complex(fun,guess,funtol=1e-3,reltol=1e-1,maxiter=10,verbose=False):
	# Expected, input a float numpy matrix guess vector of shape (L,1)
	# Expected: tol as a float
	nloop=0
	nconverg=0
	sol=n.asarray(guess)
	print sol
	nsys=len(sol)
	fx=fun(n.asarray(sol))
	x1=sol.copy()
# 	print x1
	J=n.mat(n.zeros((nsys,nsys),complex))
	normx=n.linalg.norm(sol)
	for j in xrange(nsys):
		x2=x1.copy()
		unitj=n.zeros((nsys,1),complex)
		unitj[j]=0.01*(1+1j)*reltol*normx
		unitj=n.mat(unitj)
#  			print unitj
		x2=x2+unitj
# 		print n.asarray(x2)
		fx2=fun(n.asarray(x2))
		dfdxj=(fx2-fx)/(x2[j]-x1[j])
		J[:,j]=dfdxj
	if n.linalg.det(J)==0:
		raise IOError('Determinant of Jacobian becomes zero')
	dx=-n.linalg.solve(J,fx)
	sol=sol+dx
	relerror=n.linalg.norm(dx)/n.linalg.norm(sol)
	funerror=n.linalg.norm(fx)
	if relerror < reltol and funerror<funtol:
		nloop=maxiter
###########
	if verbose:
		print "Broyden iteration on function ", fun ,"started"
	while nloop<maxiter:
		fxold=fx.copy()
		fx=fun(sol)
		df=fx-fxold
		denom=n.linalg.norm(dx)**2
		J=J+(df-J*dx)*dx.T/denom
		dx=-n.linalg.solve(J,fx)
		sol=sol+dx
		error=n.asarray(dx)
		nloop+=1
		nconverg+=1
		if nloop==maxiter-1:
			print "Warning: maximum number of iterations reached, but convergence is not reached!"
		relerror=n.linalg.norm(dx)/n.linalg.norm(sol)
		funerror=n.linalg.norm(fx)
		if verbose==True:
			print "Iteration: %g. Relative error: %.2e. Function error: %.2e" %(nloop,relerror,funerror)
# 		print "Iteration ,",nloop,", current relative error:",relerror, ", current func. value norm: ",funnorm
		if relerror<reltol and funerror<funtol:
			nloop=maxiter

	if verbose:
		print "Iteration on function ", fun ,"converged"
	return n.array(sol),nconverg,relerror
############################################################

def new_raph_list(fun,guess,tol,maxiter,verbose=False):
	nloop=0
	nconverg=0
	sol=guess.tolist()
	if type(sol)==type(1.0) or type(sol)==type(1+1j):
		nsys=1
	else:
		nsys=len(sol)
# 	print "nsys equals:",nsys
	if verbose:
		print "Newton-raphson iteration on function ", fun ,"started"
	while nloop<maxiter:
		J=n.zeros((nsys,nsys),complex)
		fx=fun(sol)
		x1=list(sol)
		for j in xrange(nsys):
			x2=list(x1)
			x2[j]=x2[j]+0.01*tol
# 			print unitj
			fx2=fun(x2)
			for i in xrange(nsys):
				dfidxj=(fx2[i]-fx[i])/(x2[j]-x1[j])
				J[i,j]=dfidxj
		print J
		if n.linalg.det(J)==0:
			raise IOError('Determinant of Jacobian becomes zero')
		dx=-n.linalg.solve(J,fx)
		sol=sol+dx
		error=n.asarray(dx)
		relerror=n.linalg.norm(dx)/n.linalg.norm(sol)
		funerror=n.linalg.norm(fx)
		if verbose:
			print "Iteration: %g. Relative error: %.3e. Function error: %.2e" %(nloop,relerror,funerror)
		nloop+=1
		nconverg+=1
		if nloop==maxiter-1:
			print "Warning: maximum number of iterations reached, but convergence is not reached!"
		abserror=n.linalg.norm((error))
		converged=False

		if abserror<tol:
			nloop=maxiter
			if verbose:
				print "Iteration on function ", fun ,"converged"
	return sol,nconverg,error


def new_raph_numpy(fun,guess,reltol,funtol=1e-6,maxiter=10,verbose=False):
	if verbose:
		print "Newton-raphson iteration on function ", fun ,"started"
	nloop=0
	nconverg=0
	sol=guess
	nsys=sol.shape[0]
	# 	print "nsys:",nsys
	normx=n.linalg.norm(guess)
# 	print "nsys equals:",nsys
	while nloop<maxiter:
		J=n.zeros((nsys,nsys),float)
		fx=fun(sol)
		x1=sol.copy()
		for j in xrange(nsys):
			x2=x1.copy()
			unitj=n.zeros((nsys,))
			unitj[j]=0.01*reltol*normx
# 			print unitj
			x2=x2+unitj
			fx2=fun(x2)
# 			print "fx:",fx
# 			print "fx2:",fx2
			dfdxj=(fx2-fx)/(x2[j]-x1[j])
# 			print "dfdxj",dfdxj
# 			print J[:,j]
			J[:,j]=dfdxj.ravel()

		if n.linalg.det(J)==0:
			print "Determinant of Jacobian equals zero!"
			k=0
			for i in xrange(nsys):
				#print j
				if sum(abs(J[i,:]))>1e-10:
					k+=1
				else:
					print "Row %g of Jacobian is (nearly) empty" %k
					k+=1


			return sol,0,0
# 			raise IOError('Determinant of Jacobian becomes zero')
		dx=(-n.linalg.solve(J,fx)).ravel()
		sol=sol+dx
		relerror=n.linalg.norm(dx)/n.linalg.norm(sol)
		funerror=n.linalg.norm(fx)
		nloop+=1
		nconverg+=1
		if verbose:
			print "Iteration: %g. Relative error: %.2e " %(nloop,relerror)
		if nloop==maxiter-1:
			print "Warning: maximum number of iterations reached, but convergence is not reached!"

		if relerror<reltol and funerror<funtol:
			nloop=maxiter
			if verbose:
				print "Iteration on function ", fun ,"converged"
				print "Relative error: %.2e. Function error: %.2e" %(relerror,funerror)
	return sol

def new_raph_num(fun,guess,reltol=1e-1,funtol=1e-3,maxiter=10,verbose=False):
	if verbose:
		print "Newton-Raphson iteration started with function tolerance %.2e and relative solution tolerance %.2e" %(reltol,funtol)
		print "Iteration on function ", fun ,"started"
	nloop=0
	nconverg=0
	sol=guess
	normx=abs(guess)
# 	print "nsys equals:",nsys
	while nloop<maxiter:
		fx=fun(sol)
		x1=sol
		x2=x1+0.01*reltol*normx
		fx2=fun(x2)
		dfdx=(fx2-fx)/(x2-x1)
		dx=-fx/dfdx
		sol=sol+dx
		nloop+=1
		nconverg+=1
		if sol!=0.0:
			relerror=abs(dx/sol)
		else:
			relerror=abs(dx)
		funerror=abs(fx)
		if verbose:
			print "Iteration: %g. Function evaluation at x=%0.3e. Relative error: %.3e. Function error: %.2e" %(nloop,sol,relerror,funerror)
		if nloop==maxiter-1:
			print "Warning: maximum number of iterations reached, but convergence is not reached!"

		if relerror<reltol and funerror < funtol:
			nloop=maxiter
			if verbose:
				print "Iteration on function ", fun ,"converged"
	return sol,nconverg,relerror

def secant(fun,guess,reltol=1e-1,funtol=1e-3,maxiter=10,verbose=False):
	if verbose:
		print "Secant iteration started with function tolerance %.2e and relative solution tolerance %.2e" %(reltol,funtol)
	nloop=0
	nconverg=0
	sol=guess
	normx=abs(guess)
# 	print "nsys equals:",nsys
	fx=fun(sol)
	x1=sol
	x2=x1+0.01*reltol*normx
	fx2=fun(x2)
	dfdx=(fx2-fx)/(x2-x1)
	dx=-fx/dfdx
	sol=sol+dx
	while nloop<maxiter:
		fxold=fx
		fx=fun(sol)
		df=fx-fxold
		denom=dx**2
		dfdx=dfdx+(df-dfdx*dx)*dx/denom
		dx=-fx/dfdx
		sol=sol+dx
		nloop+=1
		nconverg+=1
		if sol!=0.0:
			relerror=abs(dx/sol)
		else:
			relerror=abs(dx)
		funerror=abs(fx)
		if verbose:
			print "Iteration: %g. Function evaluation at x=%0.3e. Relative error: %.3e. Function error: %.2e" %(nloop,sol,relerror,funerror)
		if nloop==maxiter-1:
			print "Warning: maximum number of iterations reached, but convergence is not reached!"

		if relerror<reltol and funerror < funtol:
			nloop=maxiter
			if verbose:
				print "Iteration on function ", fun ,"converged"
	return sol,nconverg,relerror
