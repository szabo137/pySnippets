"""
module provides a wrapper for numpy evaluation of a function under a given condition
"""
import numpy as np
def condition(cond,alter=0.0):
    def condEvaluator(func):
        def tempFunc(el):
            condArr = cond(el)
            resArr=np.ones(el.shape)*alter
            resArr[condArr]=func(el[condArr])
            return resArr
        return tempFunc
    return condEvaluator
        
        
if __name__=='__main__':
    import time 
    @condition(lambda x:x>0,5.0)
    def sq(x):
        return np.sqrt(x)
    
    args = np.linspace(-5,5,20000)
    start = time.time()
    res = sq(args)
    end = time.time() - start
    print "deco: %s (time: %1.2e)"%(res,end)

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
