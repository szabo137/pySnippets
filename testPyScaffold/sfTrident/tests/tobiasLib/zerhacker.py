# -*- coding: utf-8 -*-

class zerhacker(object):
  
    def __init__(self , slicearguments = (0,) , slicelength = 50 , axis = 0):
        # this argument of self.func shall be sliced, by default it is the first
        self.sa       = slicearguments
        # size of the slices
        self.sl       = slicelength
        # relevant axis of the sliced argument
        self.axis     = axis


    def __get_slicedargument__(self,arg,slice0):
	#print '...'
	#import sys,os
	#sys.path.append( os.pardir )
	#from NME.minkowskispace import MinkowskiVector
	return arg[slice0]
	#if isinstance(arg,MinkowskiVector):

	  ##print '---> MinkowskiVector'
	  #R = MinkowskiVector( [arg[i][slice0] for i in [0,1,2,3] ] )
	  ##print R.data.shape
	  ##print R.shape
	  #return R
	#elif isinstance(arg,parray):
	  ##print '--> parray'
	  #R = arg[slice0]
	  ##print R.shape
	  #return R
	  

    def __subresult__(self , N , lower , upper , args , kwargs ):
        
	print N, lower,upper#,self.sa
        
	
        slice0            = [slice(None,None,None)] * N
        slice0[self.axis] = slice(lower,upper,None)
        slice0            = tuple(slice0)
    
	#print slice0
    
        nargs           = list(args)
	
	for sliceargument in self.sa:
	  nargs[sliceargument]  = self.__get_slicedargument__(args[sliceargument],slice0)
        nargs           = tuple(nargs)
	
	#print nargs[0].shape,nargs[1].shape
	#raw_input('__subresult__')
        
        return self.func(*nargs,**kwargs)



    def __combine_result__(self,result):
	import numpy as np

	# different combination mechanisms
	if isinstance(result[0],np.ndarray):
	    # for numpy.ndarrays
	    return self.__combine_array__(result)
	elif isinstance( result[0] , tuple ):
	    # for tuples
	    RESULT = zip(*result)
	    return tuple( map( lambda R: self.__combine_array__(R) , RESULT ) )

	

    def __combine_array__(self,result):
	import numpy as np
	if isinstance(result[0],np.ndarray):
	  # concatenate if array
	  return np.concatenate( tuple(result) , axis = self.axis )
	  # do nothing if not an array
	else:
	  return result


    def __wrapper__(self , *args , **kwargs ):
      
        # the arguments which shall be sliced
        A               = [args[ i ] for i in self.sa ]
        # the shape of the arguments which shall be sliced
        SA              = [a.shape for a in A]
        # the number of dimensions of the array A
        D               = [len( sa ) for sa in SA]
	# the types of the arguments which shall be sliced
	#T		= [ a.__class__ for a in A]
	#self.datatypes	= T
	
	#print SA
	#print T
	
        # the relevant axis of A
        axis            = self.axis
        #axis            = A.axis
	
	# assume that all sliced arguments have the same shape
        
        if SA[0][ axis ] > self.sl:
            result     = [self.__subresult__( D[0] , n*self.sl     , (n+1)*self.sl  ,args,kwargs) for n in range(SA[0][ axis ]/self.sl)]
            result.append(self.__subresult__( D[0] , (n+1)*self.sl , SA[0][ axis ]     ,args,kwargs))
            
	    #print SA
	    #print self.sl
	    ##print SA[ self.sl ]
	    #print SA[ axis ]
	    #print len(result)
	    #print [r.shape for r in tuple(result)]
	    #print tuple(result)[0]
	    
            return self.__combine_result__( tuple( result ) )
        else:
            return self.func(*args,**kwargs) 
      
    def __call__(self,func):
            
        self.func               = func
        self.__name__           = '<slicer> @ %s' % self.func.__name__
        	
        return self.__wrapper__
	
	
if __name__ == '__main__':
  import sys,os
  sys.path.append( os.pardir)
  from scipy import *
  from asimps import *
  from NME import *
  
  A	= MinkowskiVector( [linspace(0,1,100)[:,newaxis],0,0,0] )
  B	= parray(linspace(0,1,120)[newaxis,:])

  AA	= A + (0*A*B)
  BB	= B + (0*A*A)
  #print isinstance(A,MinkowskiVector)
  print AA.__class__
  print AA.shape
  print BB.shape
  raw_input('')

  f	= lambda x,A,B: (A*A)*sin(B*x)
  
  @slicer(slicearguments = (0,1) , slicelength = 7 , axis = 1)
  @slicer(slicearguments = (0,1) , slicelength = 5 , axis = 0)
  def integrate(A,B,a,b):
    return asimps(lambda x: f(x,A,B) , a, b)
  
  I= integrate(AA,BB,0,1)
  print I
  print I.shape
  