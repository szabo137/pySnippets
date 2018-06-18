# -*- coding: utf-8 -*-
from numpy import *
from parray import *

### NME - minkowskispace #################################################################
# v1.0:	  2010-12-13	Daniel Seipt:
# v1.1:   2011-01-17	Tobias Nousch:	Lichtkegelkomponenten als Methoden
# v1.1.1: 2011-02-01	Daniel Seipt:	Lichtkegelkomponenten mit kürzeren Namen, Rückgabe als parray
# v1.1.2: 2011-02-14	Daniel Seipt:	Korrektur in methode isonshell: funktioniert nun auch für m = 0
# v1.1.3: 2011-02-14	Daniel Seipt:	shape gibt nur den parray shape zurück, also ohne die erste 4
# v1.1.4: 2011-02-24	Daniel Seipt:	mass shell fest in Klasse verankert
# v2.0:   2011-03-08	Daniel Seipt:	Verhalten von getitem geändert, Komponenten des Vierervektors nun mit ._0() etc.

def isscalar(x):
    return isinstance(x,int) or isinstance(x,float) or isinstance(x,complex)


class MinkowskiVector(object):
	
	def __init__(self,A,mass = None):
		
		if (isinstance(A,ndarray) or isinstance(A,list) or isinstance(A,tuple)) and len(A) == 4:
			
			self.isparray	= any([isinstance(a,parray) for a in A])
			
			try:
				s	= (A[0]*A[1]*A[2]*A[3]).shape
				o	= ones(s)
				
				A	= [a*o for a in A]
				
			except:
				pass
			
			if not self.isparray:
				self.data	= asarray(A)
			else:
				# make shapes fit
				shapes		= [asarray(a).shape for a in A if isinstance(a,parray)]
				fitA		= [ a*ones(shapes[0]) for a in A]
				self.data	= asarray(fitA)
				
			self.shape	= self.data.shape[1:]
			self.mass	= mass
		else:
			raise TypeError	
	

	def getslices(self,othershape):
		# shape of parameter arrays (== for MinkowskiVectors without first index)
		selfshape		= self.data.shape[1:]
		#selfshape		= self.shape
		selfN,otherN		= len(selfshape),len(othershape)
		
		#print selfshape,othershape
		#raw_input('###')
		
		coincide	= lambda s,o: (s == o) or (s == 1) or (o == 1)
		
		if selfN == 0:
			selfind		= (None,)        * otherN
			otherind	= (slice(None),) * otherN

		elif otherN == 0:
			selfind		= (slice(None),) * selfN
			otherind	= (None,)        * selfN
			
			
		elif (selfN == otherN):
			if all(coincide(s,o) for s,o in zip(selfshape,othershape)):
				selfind		= (slice(None),) * selfN
				otherind	= (slice(None),) * otherN
			else:
				raise TypeError,"Parameter arrays have different shapes, cannot be multiplied uniquely"
		
		else:
			raise TypeError,"Parameter arrays have different shapes, cannot be multiplied uniquely"
		
		return selfind,otherind
		
	def __getitem__(self,index):
		return MinkowskiVector( map(lambda i: self.data[i][index] , [0,1,2,3] ) )
		#return parray(self.data.__getitem__(index))
		
	def __mul__(self,other):
		if isinstance(other,MinkowskiVector):
			# Scalar Product of 2 MinkowskiVectors

			selfindbase,otherindbase	= self.getslices(other.data.shape[1:])
			
			selfind_time	= (0,)             + selfindbase
			selfind_space	= (slice(1,None),) + selfindbase
			otherind_time	= (0,)             + otherindbase
			otherind_space	= (slice(1,None),) + otherindbase
			
			return parray(self.data[selfind_time]*other.data[otherind_time]	\
					- sum(self.data[selfind_space] * other.data[otherind_space], axis=0))
	
		elif isscalar(other):
			# Product of MinkowskiVector and Scalar
			return MinkowskiVector(other * self.data)
			
		elif isinstance(other,parray):
			#print other
			#print other.shape
			#raw_input('__mul__')
			
			selfindbase,otherindbase	= self.getslices(other.shape)
			
			selfind			= (slice(None),) + selfindbase
			otherind		= (None,)        + otherindbase
	
			return MinkowskiVector(self.data[selfind] * other[otherind])
			

		
		else:
			raise TypeError,"Error: Multiplication of %s * %s is not defined" % (self.__class__,other.__class__)

	
	def __rmul__(self,other):
		return self * other
		
	def __div__(self,other):
		if isscalar(other) or isinstance(other,parray):
			# Division of MinkowskiVector and Scalar
			return self * (1/ other)
		
		else:
			raise TypeError, "Error: Division of %s / %s is not defined" % (self.__class__,other.__class__)
	
	
	def __add__(self,other):
		if isinstance(other,MinkowskiVector):
			# sum of 2 MinkowskiVectors
			
			selfindbase,otherindbase	= self.getslices(other.data.shape[1:])
			
			selfind				= (slice(None),) + selfindbase
			otherind			= (slice(None),) + otherindbase
			
			return MinkowskiVector(self.data[selfind] + other.data[otherind])
			
		else:
			raise TypeError,"Error: Sum of %s + %s is not defined" % (self.__class__,other.__class__)
	
		
	def __neg__(self):
		return (-1)*self	
		
	def __sub__(self,other):
		if isinstance(other,MinkowskiVector):
			# difference of 2 MinkowskiVectors
			return self + (-1.*other)
		else:
			raise TypeError, "Error: Sum of %s + %s is not defined" % (self.__class__,other.__class__)

	def isonshell(self,mass = None ,tolerance = 1e-3):
		# check if MinkowskiVector is on mass shell with tolerance
		
		if mass == None and self.mass == None:
		    raise TypeError, "Error: mass shell value not given"
		elif mass == None:
		    m	= self.mass
		else:
		     m	= mass
		  
		  
		zero	= self * self - m**2
		#print massshell
		#print abs(massshell)       < tolerance
		#raw_input()
		if m == 0:	return all(abs(zero) < tolerance       )
		else:		return all(abs(zero) < tolerance * m**2)
	
	def conjugate(self):
		return MinkowskiVector(conjugate(self.data))

	def __call__(self):
		return self.data 
		
	
	def __repr__(self):
		return str(self.data)
		
	
	def component(self,i):
		return parray(self.data[i])
		
	def _0(self):
		return parray(self.data[0])
		
	def _1(self):
		return parray(self.data[1])
	
	def _2(self):
		return parray(self.data[2])
	
	def _3(self):
		return parray(self.data[3])
	    
	

	def plus_component(self):
		return parray(self.data[0]+self.data[3])

	def minus_component(self):
		return parray(self.data[0]-self.data[3])

	def perp_component(self):
		A1 = parray(self.data[1])
		A2 = parray(self.data[2])
		return [A1,A2]
		
	minus	= minus_component
	plus	= plus_component
	
	def __take__(self,index):
	  
	  return MinkowskiVector( map(lambda i: self.data[i][index] , [0,1,2,3] ) )



if __name__ == '__main__':
	
	P	= MinkowskiVector([2.354,1.112,0.1286,-0.73513])
	Q	= MinkowskiVector([0.125,1.565,-1.564,10.73513])
	print P.__class__
	print P
	print (P * P) * Q
	print P * (P * Q)
	
	print P.plus(), P._0() + P._3()
	#raw_input()
	print P.isonshell(2.)
	
	
	theta	= parray(linspace(0,pi,30))
	

	
	E0	= [1,sin(theta),0,cos(theta)]
	#print E0
	#for e in E0: print e
	#print E0.shape
	E	= MinkowskiVector(E0)
	
	print E
	print E*E
	
	print E/2
	
	
	print theta.__class__
	print isinstance(theta,ndarray)
	print theta*E
	print E*theta
	print P.plus()
	
	print '+++++++++++++++++++++++++'
	# test __take__ method
	
	print E._0()
	print E._1()
	print E[5:].shape
	#print E.shape,E.__class__
	#E2 =  E.__take__(slice(None,None,2))
	#print E2 
	#print E2.shape
	#print E2.__class__
	
	
