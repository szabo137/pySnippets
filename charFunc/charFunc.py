"""
module provides a wrapper for numpy evaluation of a function under a given condition
"""
import numpy as np
from functools import wraps
try:
    from numba import jit,vectorize,float64
    fast = jit('float64(float64[:])')
except ImportError:
    fast= lambda x:x

def condition(cond,alter=0.0):
    def condEvaluator(func):
        def tempFunc(el):
            condArr = cond(el)
            resArr=np.ones(el.shape)*alter
            resArr[condArr]=func(el[condArr])
            return resArr
        return tempFunc
    return condEvaluator

class condition2(object):
    def __init__(self,cond,alter=0.0):
        self.__cond = cond
        self.__alter=alter
    
    def __evaluator(self,func):
        @wraps(func)
        def tempFunc(el):
            condArr = self.__cond(el)
            resArr=np.ones(el.shape)*self.__alter
            resArr[condArr]=func(el[condArr])
            return resArr
        return fast(tempFunc)
    
    def __call__(self,func):
        return self.__evaluator(func)

def condition3(cond,alter=0.0):
    def condEvaluator(func):
        @wraps(func)
        def tempFunc(el):
            condArr = cond(el)
            resArr=np.ones(el.shape)*alter
            resArr[condArr]=func(el[condArr])
            return resArr
        return tempFunc
    return condEvaluator

if __name__=='__main__':
    import time 
    @condition(lambda x:x>0,0.0)
    def sq(x):
        return np.sqrt(x)*np.sqrt(x)*np.sqrt(x)*np.sqrt(x)*np.sqrt(x)
    
    args = np.linspace(-5,5,200000)
    start = time.time()
    res = sq(args)
    end = time.time() - start
    print "deco1: %s (time: %1.2e)"%(res,end)

    @condition2(lambda x:x>0,0.0)
    def sq(x):
        return np.sqrt(x)*np.sqrt(x)*np.sqrt(x)*np.sqrt(x)*np.sqrt(x)
    
    args = np.linspace(-5,5,200000)
    start = time.time()
    res = sq(args)
    end = time.time() - start
    print "deco2: %s (time: %1.2e)"%(res,end)
    
    @condition3(lambda x:x>0,0.0)
    def sq(x):
        return np.sqrt(x)*np.sqrt(x)*np.sqrt(x)*np.sqrt(x)*np.sqrt(x)
    
    args = np.linspace(-5,5,2)
    start = time.time()
    res = sq(args)
    end = time.time() - start
    print "deco3: %s (time: %1.2e)"%(res,end)

    def sqpy(x):
        tres = np.zeros(x.shape)
        for index in np.arange(len(x)):
            if x[index]>0:
                tres[index]=np.sqrt(x[index])
        return tres

    start = time.time()
    res = sqpy(args)
    end = time.time() - start
    print "py: %s (time: %1.2e)"%(res,end)

    def sq2(x):
        if x>=0:
            return np.sqrt(x)
        else:
            return 0.0
    
    
    start = time.time()
    res = np.array(map(sq2,args))
    end = time.time() - start
    print "map: %s (time: %1.2e)"%(res,end)

    def sq2(x):
        if x>=0:
            return np.sqrt(x)
        else:
            return 0.0
    sq3 = np.vectorize(sq2)
    
    start = time.time()
    res = sq3(args)
    end = time.time() - start
    print "vec: %s (time: %1.2e)"%(res,end)
    
    
    from numba import vectorize, float64
    import math as m
    
    @vectorize([float64(float64)])
    def mySq(x):
        return m.sqrt(x)
    
    @vectorize([float64(float64)])
    def sq4(x):
        if x>=0:
            return np.sqrt(x)
        else:
            return 0.0
    
    start = time.time()
    res = sq4(args)
    end = time.time() - start
    print "numba: %s (time: %1.2e)"%(res,end)
