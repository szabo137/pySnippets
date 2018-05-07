"""
test product of arrays of mks vectors
"""
"""
prodSum of mks vectors from qft lib
"""
import numpy as np
import qft
import itertools as it


def oldProdSingle(arr1,arr2):
    return np.array([arr1[item]*arr2[item2] for item in np.arange(4) for item2 in np.arange(4)]).reshape(16,)
    
def oldProd(arr1,arr2):
    return np.array([oldProdSingle(arr1[item],arr2[item]) for item in np.arange(4) for item2 in np.arange(4) ])
    
def iterProd(arr1,arr2):
    return np.prod(np.array([t for t in it.product(arr1,arr2)]),axis=1)
    
def numpyProdSingle(arr1,arr2):
    return np.multiply.outer(arr1,arr2).flatten()

def numpyProd(arr1,arr2):
    return np.array([numpyProdSingle(arr1[item],arr2[item]) for item in np.arange(4) for item2 in np.arange(4) ])
    
def numpyProdOuter(arr1,arr2):
    return np.multiply.outer(arr1,arr2)
    
    
if __name__=='__main__':
    import time as T
    testArr11 = np.array([qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([2,0,0,0]),qft.MinkowskiVector([3,0,0,0]),qft.MinkowskiVector([4,0,0,0])])
    testArr12 = np.array([qft.MinkowskiVector([2,0,0,0]),qft.MinkowskiVector([2,0,0,0]),qft.MinkowskiVector([3,0,0,0]),qft.MinkowskiVector([4,0,0,0])])
    testArr13 = np.array([qft.MinkowskiVector([3,0,0,0]),qft.MinkowskiVector([2,0,0,0]),qft.MinkowskiVector([3,0,0,0]),qft.MinkowskiVector([4,0,0,0])])
    testArr14 = np.array([qft.MinkowskiVector([4,0,0,0]),qft.MinkowskiVector([2,0,0,0]),qft.MinkowskiVector([3,0,0,0]),qft.MinkowskiVector([4,0,0,0])])
    testArr1 = np.array([testArr11,testArr12,testArr13,testArr14])
    testArr21 = np.array([qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0])])
    testArr22 = np.array([qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0])])
    testArr23 = np.array([qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0])])
    testArr24 = np.array([qft.MinkowskiVector([2,0,0,0]),qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0]),qft.MinkowskiVector([1,0,0,0])])
    testArr2 = np.array([testArr21,testArr22,testArr23,testArr24])
    
    start = T.time()
    resOld = oldProd(testArr1,testArr2)#.reshape(4,4,16)
    end = T.time() - start
    print"resOld: %s\n\ttime: %1.2e"%(resOld,end)
    print"resOld shape: %s"%(str(resOld.shape))
    
    """
    start = T.time()
    prodResParray = 25*resOld
    end=T.time() - start
    print"prodParray: %s\n\ttime: %1.2e"%(prodResParray,end)
    """
    """
    start = T.time()
    floatRes = np.asfarray(resOld)
    end=T.time() - start
    print"resFloat: %s\n\ttime: %1.2e"%(floatRes,end)
    """
    """
    start = T.time()
    prodResNP = 25*resOld
    end=T.time() - start
    print"prodResNP: %s\n\ttime: %1.2e"%(prodResNP,end)
    """
    """
    start = T.time()
    resIter = iterProd(testArr1,testArr2)
    end = T.time() - start
    print"resIter: %s\n\ttime: %1.2e"%(resIter,end)
    """
    
    start = T.time()
    resNumpy = numpyProd(testArr1,testArr2)
    end = T.time() - start
    print"resNumpy: %s\n\ttime: %1.2e"%(resNumpy,end)
    print"resNumpy shape: %s"%(str(resNumpy.shape))
    
    start = T.time()
    resNumpyOut = np.asfarray(numpyProdOuter(testArr1,testArr2))
    end = T.time() - start
    print"resNumpyOuter: %s\n\ttime: %1.2e"%(resNumpyOut,end)
    
    
    test1 = qft.MinkowskiVector([1,0,0,0])
    test2 = qft.MinkowskiVector([1,0,0,0])
    start = T.time()
    for item in np.arange(16):
        test1*test2

    end=T.time() - start
    print "single prod time: %1.2e"%(end)
