# -*- coding: utf-8 -*-
import sys,os


from parray import *
from minkowskispace import *
from numpy import *


class densitymatrix(object):
	
	
	def __init__(self,S,pshape):
		# M must be a matrix element, data type array, size (2x2) representing the outgoing photon polarizations
		# each element of M should be of data type parray (not yet tested)
		#
		# the spin summation has to be done before the normalization.
		# otherwise the spin-flip very small contribtuions in the matrix element are enhanced to the same
		# size as the non-spin-flip contributions!!!

		

		self.shape 	= (4,4) + pshape

		if isinstance(S,ndarray):
			if S.shape == (2,2,2,2):	self.data	= self.construct_density_matrix(S)
			else:				raise TypeError
		elif isinstance(S,list):
			if any(S_l.shape == (2,2,2,2) for S_l in S):	self.data	= sum( [self.construct_density_matrix(S_l) for S_l in S] , axis = 0)
			else:						raise TypeError
		else:
			raise TypeError

		self.normalize()


	def normalize(self):
		# Normierung
		rho0		= self.data
		normalization	= sum((rho0[i,i,...] for i in range(4)),axis=0)
		self.data	= rho0 / normalization[newaxis,newaxis,...]


	def construct_density_matrix(self,S):
			
			# Abbildung der 2 Polarisationsindices auf eine Achse der Dichtematrix
			ind	= self.mapping()
			S0	= asarray([asarray(S[ind[i]]) for i in range(4)])
			


			# S0.shape = (4,2,2)
			# Anzahl der Parameterdimensionen, jedes Element von S0 ist vom Typ parray
			psig	= S0[0,0,0].shape
			#psiglen	= len(psig)
			


			signature1	= (None,) + (slice(None,None,None),) + (slice(None,None,None),) * 2
			signature2	= (slice(None,None,None),) + (None,) + (slice(None,None,None),) * 2

			rho00		= sum(sum(S0[signature1]*conjugate(S0[signature2]) ,axis = -1),axis=-1)
			rho0		= zeros((4,4) + psig, dtype=complex)


			for i in range(4):
				for j in range(4):
					rho0[i,j,...]	= rho00[i,j]

			return rho0
					


	def mapping(self):
		# defines the mappings of M-indices to the mapping of the indices of the density matrix
		
		return ( ( slice(None,None,None) , 0 , 0 , slice(None,None,None) ), 
			 ( slice(None,None,None) , 0 , 1 , slice(None,None,None) ),
			 ( slice(None,None,None) , 1 , 0 , slice(None,None,None) ),
			 ( slice(None,None,None) , 1 , 1 , slice(None,None,None) )   )
	
	def trace(self):
		return sum((self.data[i,i,...] for i in range(4)),axis=0)#.real
	
	
	def sigma2(self):
		# define sigma2 (kronecker) sigma2
		return asarray([[0,0,0,-1],[0,0,1,0],[0,1,0,0],[-1,0,0,0]],dtype=float)
		

	def __getitem__(self,ind):
		return self.data[ind]
	
	def Q(self):
		
		siglen	= len(self.data.shape) - 2	# Anzahl der Dimensionen fuer das parray
		
		signature1	= (slice(None,None,None),) * 2 + (None,) * 3            + (slice(None,None,None),) * siglen
		signature2	= (None,) + (slice(None,None,None),) * 2 + (None,) * 2  + (None,) * siglen
		signature3	= (None,) * 2 + (slice(None,None,None),) * 2 + (None,)  + (slice(None,None,None),) * siglen
		signature4	= (None,) * 3 + (slice(None,None,None),) * 2            + (None,) * siglen
		
		
		sigma2	= self.sigma2()

		
		Q0	= self.data[signature1] * sigma2[signature2] * conjugate(self.data)[signature3] * sigma2[signature4]
	  
		Q	= sum(sum(sum(Q0,axis=3),axis=2),axis=1)
	
		self.Qmatrix	= Q
		return Q
	
	
	def concurrence(self):
		#import scipy
		from scipy.linalg import eig
		
		Q		= self.Q()
		parrayshape	= self.shape[2:]
		
		X		= map(range,parrayshape)
		
		if len(X) == 2:
			indices	= [(x0,x1) for x0 in X[0] for x1 in X[1]]
		elif len(X) == 3:
			indices	= [(x0,x1,x2) for x0 in X[0] for x1 in X[1] for x2 in X[2]]
		
		eigenvalue_array	= zeros((4,)+parrayshape, dtype = complex128)
		eigenvector_array	= zeros((4,4)+parrayshape, dtype = complex128)
		concurrence		= zeros(parrayshape, dtype = complex128)
		
		#print indices
		for ind in indices:
			
			eigenvalues,eigenvectors	= eig( Q[(slice(None,None,None),)*2 + ind] )

			#if ind == (0,0):
				#print ind
				#from scipy.io import write_array
				#write_array('Qmatrix_real.dat',Q[(slice(None,None,None),)*2 + ind].real,precision=16)
				#write_array('Qmatrix_imag.dat',Q[(slice(None,None,None),)*2 + ind].imag,precision=16)
	
				#write_array('Qeigenvalues_real.dat',eigenvalues.real,precision=16)
				#write_array('Qeigenvalues_imag.dat',eigenvalues.imag,precision=16)
				
				#write_array('Qeigenvectors_real.dat',eigenvectors.real,precision=16)
				#write_array('Qeigenvectors_imag.dat',eigenvectors.imag,precision=16)
				
				#A	= sum(Q[(slice(None,None,None),)*2 +(None,)+ ind] * eigenvectors[newaxis,:,:],axis=1) - eigenvalues[newaxis,:] * eigenvectors
				#print A
				
				#write_array('Qeveq_real.dat',A.real,precision=16)
				#write_array('Qeveq_imag.dat',A.imag,precision=16)	
				
				#self.q1	= Q[(slice(None,None,None),)*2 + ind]
				#self.e1	= eigenvalues
				#self.v1	= eigenvectors
				
				
			eigenvalue_array[(slice(None,None,None),) + ind]	= eigenvalues
			eigenvector_array[2*(slice(None,None,None),) + ind]	= eigenvectors

			sqrtev	= map(lambda z: sqrt(z),eigenvalues)
			zeta	= sort(sqrtev)
		
			
		
			concurrence[ind]= max(0,zeta[3] - sum(zeta[:3]))
			#print zeta,concurrence[ind]
		self.eigenvalues	= eigenvalue_array
		self.eigenvectors	= eigenvector_array
		#sys.exit()
		return concurrence
			
			
	#def concurrence_backup(self):
		#from scipy.linalg import eig
		
		#Q		= self.Q()
		#parrayshape	= self.shape[2:]
		
		#X		= map(range,parrayshape)
		
		#if len(X) == 2:
			#indices	= [(x0,x1) for x0 in X[0] for x1 in X[1]]
		#elif len(X) == 3:
			#indices	= [(x0,x1,x2) for x0 in X[0] for x1 in X[1] for x2 in X[2]]
		
		#eigenvalue_array= zeros(parrayshape+(4,), dtype = float)
		#concurrence	= zeros(parrayshape, dtype = float)
		
		#for ind in indices:
			#eigenvalues	= eig(asarray(Q[(slice(None,None,None),)*2 + ind]))[0]
			
			#sqrtev	= map(lambda z: abs(sqrt(z)),eigenvalues)
			#zeta	= sort(sqrtev)
			
			#concurrence[ind]= max(0,zeta[3] - sum(zeta[:3]))
		
		#return concurrence
	
	def entanglement(self):
		
		# definition of entanglement of formation from  Plenio arXiv:quant-ph/0504163v3 (2006)
		C	= self.concurrence()
		s	= lambda x: -x*log(x)/log(2.) - (1-x)*log(1-x)/log(2.)
		
		#x	= linspace(0+1e-6,1-1e-6,500)
		#import pylab
		#S=s( (1+sqrt(1-x**2)) / 2.)
		#pylab.plot(x,S)
		#print S
		#raw_input()
		return s( (1+sqrt(1-C**2)) / 2. )		






class densitymatrix_backup(object):
	
	
	def __init__(self,S):
		# M must be a matrix element, data type array, size (2x2) representing the outgoing photon polarizations
		# each element of M should be of data type parray (not yet tested)
		#
		#from scipy.linalg import eig
		#self.eig	= eig
		
		if not S.shape == (2,2,2,2):
			raise TypeError
		else:
			# Abbildung der 2 Polarisationsindices auf eine Achse der Dichtematrix
			ind	= self.mapping()
			S0	= asarray([asarray(S[ind[i]]) for i in range(4)])
			
			
			
			# S0.shape = (4,2,2), die letzten 2 Achsen sind die Spineinstellungen
			# Anzahl der Parameterdimensionen, jedes Element von S0 ist vom Typ parray
			siglen	= len(S0[0,0,0].shape)
			
			
			#
			signature1	= (None,) + (slice(None,None,None),) + (slice(None,None,None),) * 2
			signature2	= (slice(None,None,None),) + (None,) + (slice(None,None,None),) * 2
			#signaturen	= (None,) + (None,)                  + (slice(None,None,None),) * siglen


			# Konstruktion einer Matrix aus dem Vektor S0
			# rho0.shape = (4,4,2,2)
			rho0		= S0[signature1]*conjugate(S0[signature2])
			
			print rho0.shape
			
			# Einsortieren von rho0 in ein array der Form (4,4) + (2,2) + parray.shape
			rhoshape	= (4,4) + (2,2) + S0[0,0,0].shape			
			rho		= zeros(rhoshape, dtype=complex)
			for i in range(4):
				for j in range(4):
					for k in range(2):
						for l in range(2):
							rho[i,j,k,l,...]	= rho0[i,j,k,l]
			
			## Normierung der Spur auf 1, vor Spinsummation
			normalization	= sum((rho[i,i,...] for i in range(4)),axis=0).real * 4
			#self.data	= sum(sum(rho / normalization[newaxis,newaxis,...],axis=2),axis=2)
			self.data	= (rho[:,:,0,0,...] / normalization[newaxis,newaxis,0,0,...] + rho[:,:,1,1,...] / normalization[newaxis,newaxis,1,1,...])*2
			
			print normalization.shape
			import pylab
			pylab.figure(1)
			pylab.subplot(2,2,1)
			pylab.contourf(normalization[0,1,...],50)
			pylab.colorbar()
			pylab.subplot(2,2,2)
			pylab.contourf(normalization[1,0,...],50)
			pylab.colorbar()
			pylab.subplot(2,2,3)
			pylab.contourf(normalization[1,1,...],50)
			pylab.colorbar()
			pylab.subplot(2,2,4)
			pylab.contourf(normalization[0,0,...],50)
			pylab.colorbar()
			
			
			#Normierung der Spur auf 1, nach Spinsummation
			rhospinav	= sum(sum(rho,axis=2),axis=2)
			normalization	= sum((rhospinav[i,i,...] for i in range(4)),axis=0).real
			#self.data	= rhospinav / normalization[newaxis,newaxis,...]
			pylab.figure(2)
			pylab.contourf(normalization,50)
			pylab.colorbar()
			raw_input()
			self.shape 	= self.data.shape
			


#class densitymatrix3(object):
	
	
	#def __init__(self,S):
		## M must be a matrix element, data type array, size (2x2) representing the outgoing photon polarizations
		## each element of M should be of data type parray (not yet tested)
		##
		#from scipy.linalg import eig
		#self.eig	= eig
		
		#if not S.shape == (2,2,2,2):
			#raise TypeError
		#else:
			## Abbildung der 2 Polarisationsindices auf eine Achse der Dichtematrix
			#ind	= self.mapping()
			#S0	= asarray([asarray(S[ind[i]]) for i in range(4)])
			
			
			
			## S0.shape = (4,2,2), die letzten 2 Achsen sind die Spineinstellungen
			## Anzahl der Parameterdimensionen, jedes Element von S0 ist vom Typ parray
			#siglen	= len(S0[0,0,0].shape)
			
			
			##
			#signature1	= (None,) + (slice(None,None,None),) + (slice(None,None,None),) * 2
			#signature2	= (slice(None,None,None),) + (None,) + (slice(None,None,None),) * 2
			#signaturen	= (None,) + (None,)                  + (slice(None,None,None),) * siglen


			## Konstruktion einer Matrix aus dem Vektor S0 und Summation ueber die Spineinstellungen
			## rho0.shape = (4,4)
			#rho0		= sum(sum(S0[signature1]*conjugate(S0[signature2]),axis=-1),axis=-1)
			
			## Einsortieren von rho0 in ein array der Form (4,4) + parray.shape
			#rhoshape	= (4,4) + S0[0,0,0].shape			
			#rho		= zeros(rhoshape, dtype=complex)
			#for i in range(4):
				#for j in range(4):
					#rho[i,j,...]	= rho0[i,j]
			
			## Normierung der Spur auf 1			
			#normalization	= sum((rho[i,i] for i in range(4)),axis=0).real
			
			
			#self.data	= rho / normalization[signaturen]
			#self.shape 	= self.data.shape
					
			
			
			
			
	#def mapping(self):
		## defines the mappings of M-indices to the mapping of the indices of the density matrix
		
		#return ( ( slice(None,None,None) , 0 , 0 , slice(None,None,None) ),
	                 #( slice(None,None,None) , 0 , 1 , slice(None,None,None) ),
			 #( slice(None,None,None) , 1 , 0 , slice(None,None,None) ),
			 #( slice(None,None,None) , 1 , 1 , slice(None,None,None) )
		       #)
	
	#def trace(self):
		#return sum((self.data[i,i] for i in range(4)),axis=0).real
	
	
	#def sigma2(self):
		## define sigma2 (kronecker) sigma2
		#return asarray([[0,0,0,-1],[0,0,1,0],[0,1,0,0],[-1,0,0,0]],dtype=int)
		

	#def __getitem__(self,ind):
		#return self.data[ind]
	
	#def Q(self):
		
		#siglen	= len(self.data.shape) - 2
		
		#signature1	= (slice(None,None,None),) * 2 + (None,) * 3            + (slice(None,None,None),) * siglen
		#signature2	= (None,) + (slice(None,None,None),) * 2 + (None,) * 2  + (None,) * siglen
		#signature3	= (None,) * 2 + (slice(None,None,None),) * 2 + (None,)  + (slice(None,None,None),) * siglen
		#signature4	= (None,) * 3 + (slice(None,None,None),) * 2            + (None,) * siglen
		
		
		#sigma2	= self.sigma2()
		
		
		
		#Q0	= self.data[signature1] * sigma2[signature2] * conjugate(self.data)[signature3] * sigma2[signature4]
		
		#return sum(sum(sum(Q0,axis=3),axis=2),axis=1)
	
	
	#def concurrence(self):
		#from scipy.linalg import eig
		
		#Q		= self.Q()
		#parrayshape	= self.shape[2:]
		
		#X		= map(range,parrayshape)
		
		#if len(X) == 2:
			#indices	= [(x0,x1) for x0 in X[0] for x1 in X[1]]
		#elif len(X) == 3:
			#indices	= [(x0,x1,x2) for x0 in X[0] for x1 in X[1] for x2 in X[2]]
		
		#eigenvalue_array= zeros(parrayshape, dtype = object)
		#concurrence	= zeros(parrayshape, dtype = float)
		
		##print indices
		
		#for ind in indices:
			##print asarray(Q[(slice(None,None,None),)*2 + ind])
			##print asarray(Q[(slice(None,None,None),)*2 + ind]).shape
			##eigenvalue_array[ind]	= eig(asarray(Q[(slice(None,None,None),)*2 + ind]))[0]
			#eigenvalues	= eig(asarray(Q[(slice(None,None,None),)*2 + ind]))[0]
			#sqrtev	= map(lambda z: abs(sqrt(z)),eigenvalues)
			#zeta	= sort(sqrtev)
			##print zeta
			#concurrence[ind]= max(0,zeta[3] - sum(zeta[:3]))
			##raw_input()
		
		#return concurrence
	






#class densitymatrix2(object):
	
	
	#def __init__(self,S):
		## M must be a matrix element, data type array, size (2x2) representing the outgoing photon polarizations
		## each element of M should be of data type parray (not yet tested)
		##
		#from scipy.linalg import eig
		#self.eig	= eig
		
		#if not S.shape == (2,2):
			#raise TypeError
		#else:
			#ind	= self.mapping()
			
			#S0	= asarray([S[ind[i]] for i in range(4)])
			
			#siglen	= len(S[0,0].shape)
			
			#signature1	= (None,) + (slice(None,None,None),) + (slice(None,None,None),) * siglen
			#signature2	= (slice(None,None,None),) + (None,) + (slice(None,None,None),) * siglen
			#signaturen	= (None,) + (None,)                  + (slice(None,None,None),) * siglen

			#rho0		= S0[signature1]*conjugate(S0[signature2])
			#normalization	= sum((rho0[i,i] for i in range(4)),axis=0).real
			
			##rho		= zeros((4,4), dtype=object)
			##indices		= [(i,j) for i in range(4) for j in range(4)]
			
			##for i,j in indices: rho[i,j] = parray(rho0[i,j]/normalization)
			
			#self.data	= rho0 / normalization[signaturen]
			#self.shape 	= self.data.shape
					
			
			
			
			
	#def mapping(self):
		## defines the mappings of M-indices to the mapping of the indices of the density matrix
		
		#return ((0,0),(0,1),(1,0),(1,1))
	
	#def trace(self):
		#return sum((self.data[i,i] for i in range(4)),axis=0).real
	
	
	#def sigma2(self):
		## define sigma2 (kronecker) sigma2
		#return asarray([[0,0,0,-1],[0,0,1,0],[0,1,0,0],[-1,0,0,0]],dtype=int)
		

	#def __getitem__(self,ind):
		#return self.data[ind]
	
	#def Q(self):
		
		#siglen	= len(self.data.shape) - 2
		
		#signature1	= (slice(None,None,None),) * 2 + (None,) * 3            + (slice(None,None,None),) * siglen
		#signature2	= (None,) + (slice(None,None,None),) * 2 + (None,) * 2  + (None,) * siglen
		#signature3	= (None,) * 2 + (slice(None,None,None),) * 2 + (None,)  + (slice(None,None,None),) * siglen
		#signature4	= (None,) * 3 + (slice(None,None,None),) * 2            + (None,) * siglen
		
		
		#sigma2	= self.sigma2()
		
		
		
		#Q0	= self.data[signature1] * sigma2[signature2] * conjugate(self.data)[signature3] * sigma2[signature4]
		
		#return sum(sum(sum(Q0,axis=3),axis=2),axis=1)
	
	
	#def concurrence(self):
		#from scipy.linalg import eig
		
		#Q		= self.Q()
		#parrayshape	= self.shape[2:]
		
		#X		= map(range,parrayshape)
		
		#if len(X) == 2:
			#indices	= [(x0,x1) for x0 in X[0] for x1 in X[1]]
		#elif len(X) == 3:
			#indices	= [(x0,x1,x2) for x0 in X[0] for x1 in X[1] for x2 in X[2]]
		
		#eigenvalue_array= zeros(parrayshape, dtype = object)
		#concurrence	= zeros(parrayshape, dtype = float)
		
		##print indices
		
		#for ind in indices:
			##print asarray(Q[(slice(None,None,None),)*2 + ind])
			##print asarray(Q[(slice(None,None,None),)*2 + ind]).shape
			##eigenvalue_array[ind]	= eig(asarray(Q[(slice(None,None,None),)*2 + ind]))[0]
			#eigenvalues	= eig(asarray(Q[(slice(None,None,None),)*2 + ind]))[0]
			#sqrtev	= map(lambda z: abs(sqrt(z)),eigenvalues)
			#zeta	= sort(sqrtev)
			##print zeta
			#concurrence[ind]= max(0,zeta[3] - sum(zeta[:3]))
			##raw_input()
		
		#return concurrence
		##print eigenvalue_array.shape
		
		##print eig(asarray(Q[:,:,0,0]))
		##eigenvalues	= eig(asarray(Q[:,:,0,3]))[0]
		
		##print eigenvalues
		


