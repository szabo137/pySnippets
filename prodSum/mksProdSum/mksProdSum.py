"""
prodSum of mks vectors from qft lib
"""
import numpy as np
import qft
import itertools as it


def oldProd(arr1,arr2):
    return np.array([arr1[item]*arr2[item2] for item in np.arange(4) for item2 in np.arange(4)]).reshape(16,)
    
def iterProd(arr1,arr2):
    return np.prod(np.array([t for t in it.product(arr1,arr2)]),axis=1)
    
def numpyProd(arr1,arr2):
    return np.multiply.outer(arr1,arr2).flatten()
    
if __name__=='__main__':
    import time as T
    testArr1 = np.array([qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([2,0,0,0]),qft.MinkowskiVector([3,0,0,0]),qft.MinkowskiVector([4,0,0,0])])
    testArr2 = np.array([qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0])])
    
    
    start = T.time()
    resOld = oldProd(testArr1,testArr2)
    end = T.time() - start
    print"resOld: %s\n\ttime: %1.2e"%(resOld,end)
    
    start = T.time()
    resIter = iterProd(testArr1,testArr2)
    end = T.time() - start
    print"resIter: %s\n\ttime: %1.2e"%(resIter,end)
    
    start = T.time()
    resNumpy = numpyProd(testArr1,testArr2)
    end = T.time() - start
    print"resNumpy: %s\n\ttime: %1.2e"%(resNumpy,end)
    
    
    test1 = qft.MinkowskiVector([1,0,0,0])
    test2 = qft.MinkowskiVector([1,0,0,0])
    start = T.time()
    for item in np.arange(16):
        test1*test2

    end=T.time() - start
    print "single prod time: %1.2e"%(end)
