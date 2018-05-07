"""
test of asimps integrator
"""
import numpy as np
import qft
from qft.asimps import asimps


#np.seterr(all='log')

def integ(func,bound):
    #tempFunc = lambda x:(func(x)-func(-x))/x
    def tempFunc(x):
        #print "x: %s"%x
        res =(func(x) - func(-x))/x
        #print "res: %s"%res
        return res
    return asimps(tempFunc, 1e-8, bound, Nx=421 , errorabs = 1e-4, maxrecur = 100)



if __name__=='__main__':
    import time

    def func(x):
        #return np.array([np.sin(x),np.sin(x)])
        return np.ones(2)*np.sin(x**3)

    """
    test = np.linspace(1,100,1000)
    (func(test) - func(-test))/test

    """
    anaRes=np.pi/3.0

    start=time.time()
    res = integ(func,200)
    end=time.time() - start
    print "res: %s (anaErr: %s time: %1.2e)"%(res,np.abs(anaRes-res[0])/anaRes,end)
