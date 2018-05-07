"""
tests the speed of mks vectors

"""
import numpy as np
import qft
import itertools as it


def mksProd(mks1,mks2):
    return mks1*mks2
    
def elProd(arr1,arr2):
    #single
    return arr1[0]*arr2[0] - arr1[1]*arr2[1] - arr1[2]*arr2[2] - arr1[3]*arr2[3]


def elProdArrFor(arr1,arr2):
    el,dim=arr1.shape
    res=np.zeros(dim)
    for item in np.arange(dim):
        res[item] = arr1[0][item]*arr2[0][item] - arr1[1][item]*arr2[1][item] - arr1[2][item]*arr2[2][item] - arr1[3][item]*arr2[3][item]
    return res

if __name__=='__main__':
    import time as T
    
    
    rounds = 100
    print ""
    print "-"*15 + "Single Tests (rounds: %s)"%rounds + "-"*15
    print ""
    
    testArr1 = np.array([1.0,0,0,0])
    testArr2 = np.array([1.0,0,0,0])
    
    testMks1 = qft.MinkowskiVector(testArr1)
    testMks2 = qft.MinkowskiVector(testArr2)
    
    
    start = T.time()
    for item in np.arange(rounds):
        mksRes=mksProd(testMks1,testMks2)
    end = T.time() - start
    print"resMks: %s\n\ttime: %1.2e"%(mksRes,end)
    
    
    start = T.time()
    for item in np.arange(rounds):
        elRes=elProd(testArr1,testArr2)
    end = T.time() - start
    print"resEl: %s\n\ttime: %1.2e"%(elRes,end)


    rounds = 1
    dim=1000000
    print ""
    print "-"*15 + "Array Tests (dim: %s rounds: %s)"%(dim,rounds) + "-"*15
    print ""
    
    baseArray1 = np.ones(dim)
    baseArray0 = np.zeros(dim)
    baseArrayNum = np.arange(dim)
    
    testArray1 = np.array([baseArray1,baseArray0,baseArray0,baseArrayNum])
    testArray2 = np.array([baseArray1,baseArray0,baseArray0,baseArrayNum])
    
    
    start = T.time()
    testMks1 = qft.MinkowskiVector(testArray1)
    testMks2 = qft.MinkowskiVector(testArray2)
    
    for item in np.arange(rounds):
        mksRes=mksProd(testMks1,testMks2)
    end = T.time() - start
    print"resMks: %s\n\ttime: %1.2e"%(mksRes,end)


    start = T.time()
    for item in np.arange(rounds):
        elRes=elProdArrFor(testArray1,testArray2)
    end = T.time() - start
    print"resEl: %s\n\ttime: %1.2e"%(elRes,end)
    
    eps=1e-10
    for item in np.arange(dim):
        if np.abs(mksRes[item] - elRes[item])> eps:
            print "fail at <%s>"%item
            break
