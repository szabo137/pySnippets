# -*- coding: utf-8 -*-

def asimps( f , xmin , xmax , Nx = 200, errorabs = 1e-2 , maxrecur = 10 , plotting = False ):
    
    """
    
    Adaptive Simpson Integration v2
    
    f        :  f(x), callable function, only one argument, should return a scalar or an numpy array
    xmin     :  lower bound of integration interval
    xmax     :  upper bound of integration interval
    Nx       :  initial number of nodes
    errorabs :  maximal allowed absolute error
    maxrecur :  maximum recursiondepth
    plotting : plot nodes in figure
    
    return ( result , errorestimate , N )
    
    result        (array/scalar): numerical approximation of integral
    errorestimate (array/scalar): estimated absolute error
    N             (integer)     : number of subintervals where maximum recursiondepth has been reached,
                                  if N > 0, the required accuracy might not be reached

    v2.0	2011-03-03	Daniel Seipt: Code rewritten, array and scalar code combined, 
                                              better error estimation taken from "Handbook of Comp. Meth. for Integration",
                                              returns error estimate now and number of failed subintervals
                                              
    
    """
    from scipy import isscalar
    
    
    try:
        T     = f( 0.0 )
    except:
        from numpy.random import uniform
        T     = f( uniform(-1,1) )
	
    try:
	fshape = T.shape
    except:
        if isscalar(T):	fshape = tuple([])
	else:		raise TypeError, 'f returns neither scalar nor array'
    
    
    return _asimps_initialize( f , xmin , xmax , Nx , errorabs , maxrecur , fshape , plotting )
    
  

def _asimps_initialize( f , xmin , xmax , Nx , errorabs , maxrecur , fshape , plotting = False):
    # entry point of adaptive integration algorithm 
    # do first integration for error estimation and start recursion
    from numpy import linspace

    # specify the x-axis 
    axis  = 0					# OK

    # make number of integration nodes odd, so that even the halves are odd
    N     = 4 * ( int(Nx) // 4 ) + 1		# OK
    
    allnew        = (None,) * (len(fshape) + 1 )
    slice0        = _tupleset( allnew , axis ,slice(None,None,None))
    x             = linspace( xmin , xmax , N )[slice0]
    
    # evaluate function at nodes
    F     = f(x)


    # define the two half intervals and perform integration on these half intervals
    all2           = (slice(None),) * (len(fshape) + 1 )
    slice1        = _tupleset( all2 , axis ,slice(None,N/2+1,None))
    slice2        = _tupleset( all2 , axis ,slice(N/2 ,None,None ))
    
    
    F1,F2 = F[slice1]                             , F[slice2]
    x1,x2 = x[slice1]                             , x[slice2]
    B1,B2 = _basic_simps( F1 , x1 , axis ) , _basic_simps( F2 , x2 , axis )
    
    # plotting the points
    if plotting:
        from numpy import ravel,reshape,prod
        import pylab
        pylab.ion()
        
        pylab.figure(plotting)
        pylab.clf()
        
        sig = (F1.shape[0],prod(F1.shape[1:]))
        pylab.plot( x1.ravel() , F1.reshape(sig) , 'r.')
        pylab.plot( x2.ravel() , F2.reshape(sig) , 'b.')
        raw_input()
    

    # call recursion routine on subintervals, pass result for error estimation
    if maxrecur == 0:
      return ( B1 + B2 , B1 + B2 , 1 )
      #return B1 + B2
    else:
      R1 = _asimps_recursion(f , F1 , x1 , B1 , errorabs/2 , maxrecur , recursiondepth = 1 , axis = axis , plotting = plotting )
      R2 = _asimps_recursion(f , F2 , x2 , B2 , errorabs/2 , maxrecur , recursiondepth = 1 , axis = axis , plotting = plotting )
      return (R1[0] + R2[0] , R1[1] + R2[1] , R1[2] + R2[2] )
      #return _recursion_vector(f , F1 , x1 , B1 , errorabs/2 , maxrecur , recursiondepth = 1 , axis = axis , plotting = plotting ) \
	   #+ _recursion_vector(f , F2 , x2 , B2 , errorabs/2 , maxrecur , recursiondepth = 1 , axis = axis , plotting = plotting )

  
  
def _asimps_recursion( f , F , x , result , errorabs , maxrecur , recursiondepth , axis , plotting = False ):
    # recursion routine for adaptive simpson integration
    # calls itself on subintervals until desired accuracy (errorabs) is achieved or recursiondepth reaches maxrecur
    from numpy import argsort,concatenate,ravel,zeros,insert,empty
    #print recursiondepth,maxrecur,'_recursion_vector',':',errorabs 
    
    # calculate positions of new nodes and evaluate f at new nodes
    dx    = (x[1:] - x[:-1])/2
    x0    = x[:-1] + dx
    F0    = f(x0)
    
  
    
    ## combine old and new nodes and function values, sort them
    #x00   = concatenate((x,x0))
    #ind   = argsort(x00 , axis = axis)
    
    
    ## reduce dimensions of x and F, they blow up upon sorting, no idea why. this is a workaround
    #r     = len(F0.shape) - 1
    #slice0 = (slice(None),) + (0,)*r + (slice(None),)*r

    #F     = concatenate((F,F0))[ind][slice0]
    #x     = x00[ind][slice0]
    
    
    newFshape = (2*F.shape[0]-1,) + F.shape[1:]
    newxshape = (2*x.shape[0]-1,) + x.shape[1:]
    
    newF = empty(newFshape , dtype = F.dtype)
    newx = zeros(newxshape , dtype = x.dtype)
    
    newF[0::2,...] 	= F
    newF[1::2,...]	= F0
    
    newx[0::2,...] 	= x
    newx[1::2,...]	= x0
    
    x	= newx
    F	= newF
  
        
    
    
    # define the two half intervals and perform integration on these half intervals
    N      = x.shape[axis]
    all2   = (slice(None),) * len(x.shape) 
    slice1 = _tupleset( all2 , axis ,slice(None,N/2+1,None))
    slice2 = _tupleset( all2 , axis ,slice(N/2 ,None,None ))

    F1,F2 = F[slice1]                             , F[slice2]
    x1,x2 = x[slice1]                             , x[slice2]
    B1,B2 = _basic_simps( F1 , x1 , axis ) , _basic_simps( F2 , x2 , axis )

    # plotting the points
    if plotting:
        from numpy import ravel,reshape,prod
        import pylab
        pylab.ion()
        
        sig = (F1.shape[0],prod(F1.shape[1:]))
        pylab.plot( x1.ravel() , F1.reshape(sig) , 'r.')
        pylab.plot( x2.ravel() , F2.reshape(sig) , 'b.')
        raw_input()
    #print any((abs(B1+B2-result)/abs(result+1e-42)).ravel() < errorabs),all((abs(B1+B2-result)/abs(result+1e-42)).ravel() < errorabs)
    #if all((abs(B1+B2-result)/abs(result+1e-42)).ravel() < errorabs):
    
    # Fehlerschätzung aus:
    # Handbook of Computational Methods for Integration
    # Prem K. Kythe
    # Michael R Schäferkotter
    # Chapman & Hall/CRC, Boca Raton (2005)
    # Kapitel 2.10, S. 91 
    errorestimate	= abs(B1+B2-result) / 15.
    
    
    if all(errorestimate.ravel() < errorabs):
        return (B1 + B2 , errorestimate , 0 )
    # or maximum recursiondepth, ascend one level of recursion
    elif recursiondepth >= maxrecur:
        print 'Warning: maximum recursiondepth reached, desired accuracy might not be achieved'
        return (B1 + B2 , errorestimate , 1)
    # or start next recursion, descend one level of recursion
    else:
      R1 = _asimps_recursion(f , F1 , x1 , B1 , errorabs/2 , maxrecur , recursiondepth + 1 , axis = axis , plotting = plotting )
      R2 = _asimps_recursion(f , F2 , x2 , B2 , errorabs/2 , maxrecur , recursiondepth + 1 , axis = axis , plotting = plotting )
      return (R1[0] + R2[0] , R1[1] + R2[1] , R1[2] + R2[2] )



def _basic_simps( y , x , axis ):
    # these parts have been taken from /scipy/integrate/quadrature.py
    
    import numpy as np
    
    #raw_input('----')
    
    # determine number of nodes
    nd  = y.shape
    N   = nd[axis]
        
    # only equidistant nodes are allowed
    dx = x[1] - x[0]
    
    # define slices for the simpson formula
    start,stop,step     = 0 , N - 2 , 2
    
    all   = ( slice(None) , ) * len(nd)
    
    slice0 = _tupleset( all , axis , slice(start  , stop  , step) )
    slice1 = _tupleset( all , axis , slice(start+1, stop+1, step) )
    slice2 = _tupleset( all , axis , slice(start+2, stop+2, step) )

    # apply simpson formula
    return dx/3.0* np.sum(  (y[slice0]+4*y[slice1]+y[slice2]) , axis = axis)
    

def _tupleset(t, i, value):
    l = list(t)
    l[i] = value
    return tuple(l)



if __name__ == '__main__':
  
  import sys,os
  sys.path.append(os.path.expanduser('~') + '/bin/python/' )

  from gadgets import benchmark
  import pylab
  pylab.ion()

  import numpy as np

  
  
  from scipy import arange,linspace,newaxis
  from scipy.integrate import quad
  from scipy.special import jv,sin
  import pylab
  
  
  @benchmark
  def aquad(f,*args,**kwargs):
    realpart    = quad(lambda x: f(x).real , *args , **kwargs )[0]
    imagpart    = quad(lambda x: f(x).imag , *args , **kwargs )[0]
    return realpart + 1j*imagpart
    
    print _tupleset(all, axis, slice(start, stop, step))

  
  #f,l,u     = lambda x: exp(1j*linspace(1,2,10)/x),0.05,5
  #f,l,u     = lambda x: jv(arange(4),x),0,50
  asimps2 = benchmark(asimps)
  
  #f,l,u     = lambda x: linspace(1,2,11)[newaxis,:]*sin(linspace(1,2,20)[:,newaxis]/x),0.001,50
  A	= linspace(1,2,20)
  f,l,u     = lambda x: sin(A/x),0.001,50
  #f     = sin
  #f     = lambda x: x
  
  #for j in range(100):
  #C             = aquad( f  , 0.005 , 5)
  #oldadd	= True
  result1       = asimps2(f  , l , u , 200 , plotting = 0 , errorabs = 1e-4 , maxrecur = 50)
  
  pylab.semilogy()
  pylab.plot( A , result1[0] , 'r-')
  pylab.plot( A , abs(result1[1]) , 'b:')
  pylab.show()
  #oldadd	= False
  #result2       = asimps2(f  , l , u , 200 , plotting = 0 , errorabs = 1e-4 , maxrecur = 50) 
  
  
  print result1[0]# -  result1[0]
  #print aquad.__analysis__()
  #print asimps.__analysis__()
  
  #numberB       = sum(number)
  
  #print 'B',resultB
  
  #print asimps2( lambda x: 1./(1-0.998*x**2) , 0 , 1 , plotting = 0 , errorabs = 1e-8 , maxrecur = 10)
  #print '(3.803756514651015)'
  
  