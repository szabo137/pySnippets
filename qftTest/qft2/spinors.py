# -*- coding: utf-8 -*-
import sys,os


from parray import *
from mks import *
from numpy import *

### NME - spinors #################################################################
# v1.0:	  2010-12-13	Daniel Seipt:
# v1.1:   2011-01-14	Tobias Nousch: 	V-Spinor implementiert
# v2.0:   2011-03-08	Daniel Seipt:	Anpassung auf neues minkowskispace
# v2.1:   2011-06-22  Tobias Nousch: helicity-basis fuer V-Spinor



def pardot(left,right):
	    # general multiplication
	    leftshape,rightshape	= left.shape,right.shape
	    leftN,rightN		= len(leftshape),len(rightshape)

	    if leftN == rightN:
		# Multiplication of 2 DiracMatrices
		leftind	= (slice(None),)*2 + (None,) + (slice(None),)*(leftN-2)
		rightind	= (None,) + (slice(None),)*rightN
		return sum(left[leftind]  *  right[rightind], axis = 1 )

	    elif leftN == (rightN+1):
		# Multiplication DiracMatrix * BiSpinor
		#leftind	= (slice(None),)*2 + (None,) + (slice(None),)*(leftN-2)
		rightind	= (None,) + (slice(None),)*rightN
		return sum( left * right[rightind] , axis = 1 )

	    elif (leftN+1) == rightN:
		# Multiplication BiSpinorBar * DiracMatrix
		leftind	= (slice(None),) + (None,) + (slice(None),)*(leftN-1)
		#rightind	= (None,) + (slice(None),)*rightN
		return sum( left[leftind] * right , axis = 0 )


def feyndagg(P):
	gamma=GammaMatrix()
	if isinstance(P,MinkowskiVector):
		M		= gamma[0]*parray(P._0()) - ( gamma[1]*parray(P._1()) + gamma[2]*parray(P._2()) + gamma[3]*parray(P._3()) )
		return M
	else:
		raise TypeError,'No MinkowskiVector'



class BiSpinor(object):

    def __init__(self,psi):
	    self.data	= psi
	    try:	self.shape = self.data.shape
	    except:	self.shape = None

    def getslices(self,othershape):

	selfshape	= self.data.shape[1:]
	selfN,otherN	= len(selfshape),len(othershape)

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


    def __add__(self,other):
        if isinstance(other,BiSpinor):
	    selfindbase,otherindbase	= self.getslices(other.data.shape[1:])

	    selfind	= (slice(None,None,None),) + selfindbase
	    otherind	= (slice(None,None,None),) + otherindbase
	    return BiSpinor(self.data[selfind] + other.data[otherind])
	else:
	    raise TypeError,"Error: Addition of %s + %s is not defined" % (self.__class__,other.__class__)


    def __sub__(self,other):
	if isinstance(other,BiSpinor):

	    selfindbase,otherindbase	= self.getslices(other.data.shape[1:])

	    selfind	= (slice(None,None,None),) + selfindbase
	    otherind	= (slice(None,None,None),) + otherindbase


	    return BiSpinor(self.data[selfindbase] - other.data[otherind])

	else:
	    raise TypeError,"Error: Subtraction of %s - %s is not defined" % (self.__class__,other.__class__)
    def __mul__(self,other):

	if isinstance(other,BiSpinorBar):
	    # Multiplication SpinorU * SpinorUBar -> DiracMatrix == something like dyadic product
	    # parameters not yet supported
	    return DiracMatrix(self.data[:,newaxis]*other.data)

	else:
	    raise TypeError,"Error: Multiplication of %s * %s is not defined" % (self.__class__,other.__class__)

    def __div__(self,other):

	if isinstance(other,parray):

	    selfindbase,otherindbase	= self.getslices(other.shape)

	    print selfindbase,otherindbase

	    selfind	= (slice(None,None,None),) + selfindbase
	    otherind	= (None,)		   + otherindbase


	    return BiSpinor(self.data[selfind] / other[otherind])

	else:
	    raise TypeError,"Error: Division of %s / %s is not defined" % (self.__class__,other.__class__)


    def __rmul__(self,other):
	    pass

    def __repr__(self):
        return str(self.data)

    def __norm__(self):
        return sum(abs(self.data))

    def __getitem__(self,index):
	    return self.data[index]


class BiSpinorBar(BiSpinor):

    def __init__(self,psi):
	BiSpinor.__init__(self,None)
        self.data	= psi
	self.shape = self.data.shape

    def getslices(self,othershape):

	selfshape	= self.data.shape[1:]
	selfN,otherN	= len(selfshape),len(othershape)

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
	    return self.data[index]


    def __mul__(self,other):
	if isinstance(other,BiSpinor):
	    # Multiplication of BiSpinorBar * BiSpinor
	    # hier darf man pardot nicht verwenden, da anhand der Anzahl der Dimensionen
	    # der arrays keine eindeutige unterscheidung
	    # zu MatrixMatrix * DiracMatrix moeglich ist

	    selfindbase,otherindbase	= self.getslices(other.data.shape[1:])

	    selfind	= (slice(None),)	+ selfindbase
	    otherind	= (slice(None),)	+ otherindbase

	    return parray(sum(self.data[selfind] * other.data[otherind],axis=0))


	elif isinstance(other,DiracMatrix):
	    # Multiplikation von BiSpinorBar * DiracMatrix

	    selfindbase,otherindbase	= self.getslices(other.data.shape[2:])

	    selfind	= (slice(None),)	+ selfindbase
	    otherind	= (slice(None),) * 2	+ otherindbase


	    return BiSpinorBar( pardot(self.data[selfind],other.data[otherind]) )

	elif isscalar(other):
	    return BiSpinorBar(other * self.data)

	else:
		raise TypeError,"Error: Multiplication of %s * %s is not defined" % (self.__class__,other.__class__)
		return None



class DiracMatrix(object):


	def __init__(self,M):


	    if isinstance(M,DiracMatrix):
	        self.data = M.data
	    elif isinstance(M,ndarray) and M.shape[:2] == (4,4):
		self.data = array(M)
	    else:
		raise TypeError,'unsupported type for DiracMatrix'
	    self.shape = self.data.shape

	def getslices(self,othershape):

		selfshape	= self.data.shape[2:]
		selfN,otherN	= len(selfshape),len(othershape)

		coincide	= lambda s,o: (s == o) or (s == 1) or (o == 1)

		if selfN == 0:
			selfind		= (None,)        * otherN
			otherind	= (slice(None),) * otherN

		elif otherN == 0:
			selfind		= (slice(None),) * selfN
			otherind	= (None,)        * selfN

		elif (selfN == otherN):
			if all(coincide(s,o) for s,o in zip(selfshape,othershape)):
				#print 'DM:',[all(coincide(s,o) for s,o in zip(selfshape,othershape))]
				selfind		= (slice(None),) * selfN
				otherind	= (slice(None),) * otherN
			else:
				raise TypeError,"Parameter arrays have different shapes, cannot be multiplied uniquely"

		else:
			print selfshape,othershape
			raise TypeError,"Parameter arrays have different shapes, cannot be multiplied uniquely"

		return selfind,otherind


	def __neg__(self):
#		return DiracMatrix(-self.data)
		return -1*self



	def __add__(self,other):

	    if isscalar(other):
	        unity	= diag([1]*4)
		selfN	= len(self.data.shape)

		index	= (slice(None),)*2 + (None,)*(selfN-2)

		return DiracMatrix(self.data + other*unity[index])

	    elif isinstance(other,parray):

		selfindbase,otherindbase	= self.getslices(other.shape)
		unity	= diag([1]*4)

		unityind	= (slice(None),) * 2  + (None,) * len(otherindbase)
		otherind	= (None,)        * 2  + otherindbase
		selfind		= (slice(None),) * 2  + selfindbase

		return DiracMatrix(self.data[selfind] + other[otherind]*unity[unityind])



	    elif isinstance(other,DiracMatrix):
		# Addition von zwei DiracMatrizen

	    	selfindbase,otherindbase	= self.getslices(other.shape[2:])

		selfind		= (slice(None),)*2 + selfindbase
		otherind	= (slice(None),)*2 + otherindbase

		return DiracMatrix(self.data[selfind] + other.data[otherind])


	    else:
		raise TypeError,"Error: Addition of %s * %s is not defined" % (self.__class__,other.__class__)


	def __radd__(self,other):
    	    return self + other

	def __sub__(self,other):
	    return self + (-1.0*other)

	def __rsub__(self,other):
	    return -1.0*self + other

	def __mul__(self,other):

	    if isscalar(other):
		# Multiplikation einer DiracMatrix mit einem Skalar
		return DiracMatrix(self.data * other)

	    if isinstance(other,parray):

	    	selfindbase,otherindbase	= self.getslices(other.shape)


		selfind		= (slice(None),)*2 + selfindbase
		otherind	= (None,)*2	   + otherindbase

		return DiracMatrix(self.data[selfind] * other[otherind])


	    elif isinstance(other,DiracMatrix):
		# Multiplikation von zwei DiracMatrizen

	    	selfindbase,otherindbase	= self.getslices(other.shape[2:])

		selfind		= (slice(None),)*2 + selfindbase
		otherind	= (slice(None),)*2 + otherindbase

		return DiracMatrix( pardot(self.data[selfind],other.data[otherind]) )


	    elif isinstance(other,BiSpinor):
		# Multiplikation von DiracMatrix mir BiSpinor

		selfindbase,otherindbase	= self.getslices(other.shape[1:])


		selfind		= (slice(None),)*2 + selfindbase
		otherind	= (slice(None),)*1 + otherindbase


		return BiSpinor( pardot(self.data[selfind],other.data[otherind]))


	    else:
		raise TypeError,"Error: Multiplication of %s * %s is not defined" % (self.__class__,other.__class__)


	def __rmul__(self,other):

	    if isscalar(other):
	    	# Multiplikation einer DiracMatrix mit einem Skalar
		return DiracMatrix(other * self.data)

	    if isinstance(other,DiracMatrix):
		# Multiplikation von zwei DiracMatrizen
		return DiracMatrix(dot(other.data,self.data))

	    elif isinstance(other,parray):

		selfindbase,otherindbase	= self.getslices(other.shape)

		selfind		= (slice(None),)*2 + selfindbase
		otherind	= (     (None),)*2 + otherindbase

		return DiracMatrix(other[otherind] * self.data[selfind])


	    else:
		raise TypeError,"Error: Multiplication of %s * %s is not defined" % (other.__class__,self.__class__)

	def __div__(self,other):
		return self * (1./other)

	def __rdiv__(self,other):
		raise TypeError,"Error: Division by DiracMatrix is not defined"

	def trace(self):
		return self.data[0,0] + self.data[1,1] + self.data[2,2] + self.data[3,3]


	def __repr__(self):
	    return str(self.data)

	def __norm__(self):
	    return sum(abs(self.data))


	def __call__(self):
	    return self.data

	def __getitem__(self,index):
	    #print index
	    return self.data.__getitem__(index)







def SpinorU((p,m),s , eigenspinor = 'sigmaz'):
    from scipy import sign

    if not isinstance(p,MinkowskiVector):
        print p
        print type(p)
        p	= MinkowskiVector(p)
    else:
        pass


    if not p.isonshell(m):
        raise TypeError,'Momentum is not on shell: p^2 = %s'%(p*p)

    if eigenspinor == 'sigmaz':
        restframe	= zeros(4)

        if   s == 1:	restframe[0] = 1;
        elif s == 2:	restframe[1] = 1;
        else:		print 'error, s != 1 or 2'

    elif eigenspinor == 'helicity':

        theta = arctan2(sqrt(p._1()**2+p._2()**2),p._3())
        phi = arctan2(p._2(),p._1())

        try:	  rshape = (4,) + (theta*phi).shape
        except: rshape = 4

        restframe		= zeros( rshape   , dtype = complex)


        if   s == 1:
            restframe[0] =  exp(-1j*phi/2) * cos(theta/2)
            restframe[1] =  exp( 1j*phi/2) * sin(theta/2)
        elif s == 2:
            restframe[0] = -  exp(-1j*phi/2) * sin(theta/2)
            restframe[1] =    exp( 1j*phi/2) * cos(theta/2)
        else:		print 'error, s != 1 or 2'
    else:		print 'error, unknown eigenspinor type'


    booster		= (feyndagg(p) + m)/sqrt((m+parray(abs(p._0()) )))

    return booster * BiSpinor(restframe)


def SpinorUBar((p,m),s , eigenspinor = 'sigmaz' ):

	u		= SpinorU((p,m),s , eigenspinor ).data
	gamma0		= GammaMatrix()[0]
	return BiSpinorBar(conjugate(u)) * gamma0


def SpinorV((p,m),s , eigenspinor = 'sigmaz'):
    from scipy import sign

    if not isinstance(p,MinkowskiVector):
        p = MinkowskiVector(p)
    else:
        pass


    if not p.isonshell(m):
        raise TypeError,'Momentum is not on shell: p^2 = %s'%(p*p)

    if eigenspinor == 'sigmaz':
        restframe	= zeros(4)

        if   s == 1:	restframe[2] = 1;
        elif s == 2:	restframe[3] = 1;
        else:  print 'error, s != 1 or 2'

    elif eigenspinor == 'helicity':
        #print 'noch nicht fertig'
        theta	= arctan2(sqrt(p._1()**2+p._2()**2),p._3())
        phi	= arctan2(p._2(),p._1())
        try:	  rshape = (4,) + (theta*phi).shape
        except: rshape = 4

        restframe		= zeros( rshape   , dtype = complex)


        if   s == 1:
            restframe[2] =  exp(-1j*phi/2) * cos(theta/2)
            restframe[3] =  exp( 1j*phi/2) * sin(theta/2)
        elif s == 2:
            restframe[2] = -  exp(-1j*phi/2) * sin(theta/2)
            restframe[3] =    exp( 1j*phi/2) * cos(theta/2)

        else:   print 'error, s != 1 or 2'
    else: print 'error, unknown eigenspinor type'


    booster		= (-feyndagg(p) + m)/sqrt(m+parray(abs(p._0()) ))

    return booster * BiSpinor(restframe)



def SpinorVBar((p,m),s , eigenspinor = 'sigmaz' ):

	v		= SpinorV((p,m),s , eigenspinor ).data
	gamma0		= GammaMatrix()[0]
	return BiSpinorBar(conjugate(v)) * gamma0



class GammaMatrix(object):

    def __init__(self):

	gamma0	= array([ [ 1, 0, 0, 0],
                    [ 0, 1, 0, 0],
                    [ 0, 0,-1, 0],
                    [ 0, 0, 0,-1]])

	gamma1	= array([ [ 0, 0, 0,-1],
                    [ 0, 0,-1, 0],
                    [ 0, 1, 0, 0],
                    [ 1, 0, 0, 0]])

	gamma2	= array([ [  0, 0 , 0 , 1j],
                    [  0, 0 ,-1j,  0],
                    [  0,-1j, 0 ,  0],
                    [ 1j, 0 , 0 ,  0]])

	gamma3	= array([ [ 0, 0,-1, 0],
                    [ 0, 0, 0, 1],
                    [ 1, 0, 0, 0],
                    [ 0,-1, 0, 0]])

	self.gammavec	= np.array([DiracMatrix(gamma0),DiracMatrix(gamma1),DiracMatrix(gamma2),DiracMatrix(gamma3)])

    def __call__(self):
	return self.gammavec

    def __getitem__(self,index):
	return self.gammavec[index]






class gamma5(DiracMatrix):

    def __init__(self):

      gamma	= GammaMatrix()()
      self.data	= (-1j*gamma[0] * gamma[1] * gamma[2] * gamma[3]).data
      self.shape = self.data.shape



def spinaverage(M):
	from numpy import sum
	return sum(sum(M,axis=-1),axis=0) / 2.
