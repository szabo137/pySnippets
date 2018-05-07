"""
module to test sums of products of 1D-arrays
"""
import numpy as np


def singleProdSumFor(arrayA,arrayB):
    res = 0.0
    for index in np.arange(len(arrayA)):
        res+=arrayA[index]*arrayB[index]
    return res

def singleProdSumNumpy(arrayA,arrayB):
    return np.sum(arrayA*arrayB)


def singleProdSumEinsum(arrayA,arrayB):
    return np.einsum("j,j",arrayA,arrayB)

def prodSumFor(arrayA,arrayB,arrayC,arrayD):
    res = 0.0
    for index in np.arange(len(arrayA)):
        res+=arrayA[index]*arrayB[index]*arrayC[index]*arrayD[index]
    return res

def prodSumNumpy(*arrays):
    """
    calculation of the product of several 1-D arrays elementwise and summation of the products
    """
    return np.sum(np.prod(arrays,axis=0))

def prodSumEinsum(*arrays):
    temp = np.array(arrays)
    mod = "".join(["j," for el in np.arange(len(temp))])[:-1]
    return np.einsum(mod,*temp)

if __name__=='__main__':
    from time import time as T
    
    def anaResCalc(num,fac):
        return (endPoint-1)*endPoint*(2*(endPoint-1) + 1)/6.0*fac
    
    endPoint = 1000
    a=np.arange(endPoint)
    b=2.0*np.arange(endPoint)
    anaRes = anaResCalc(endPoint,2)
    
    print"============ 4 arrays (dim=%s) ================"%endPoint
    print "analytical res: %s"%anaRes   
    
    #for loop
    start = T()
    resFor = singleProdSumFor(a,b)
    end = T() - start
    print"for loop: %s (time: %1.2e)"%(resFor,end)
    
    #numpy
    start = T()
    resNumpy = singleProdSumNumpy(a,b)
    end = T() - start
    print"numpy: %s (time: %1.2e)"%(resNumpy,end)
    
    #einsum
    start = T()
    resEinsum = singleProdSumEinsum(a,b)
    end = T() - start
    print"einsum: %s (time: %1.2e)"%(resEinsum,end)

    A=np.arange(endPoint)
    B=np.arange(endPoint)
    C=2*np.ones(endPoint)
    D=2*np.ones(endPoint)
    
    print"\n============ 4 arrays (dim=%s) ================"%endPoint
    
    anaRes = anaResCalc(endPoint,4)
    print "anaRes: %s"%(anaRes)
    
    #for loop
    start = T()
    resFor = prodSumFor(A,B,C,D)
    end = T() - start
    print"for loop: %s (time: %1.2e)"%(resFor,end)
    
    #numpy
    start = T()
    resNumpy = prodSumNumpy(A,B,C,D)
    end = T() - start
    print"numpy: %s (time: %1.2e)"%(resNumpy,end)
    
    #einsum
    start = T()
    resEinsum = prodSumEinsum(A,B,C,D)
    end = T() - start
    print"einsum: %s (time: %1.2e)"%(resEinsum,end)

    A=np.arange(endPoint)
    B=np.arange(endPoint)
    C1=2*np.ones(endPoint)
    C2=2*np.ones(endPoint)
    C3=2*np.ones(endPoint)
    C4=2*np.ones(endPoint)
    C5=2*np.ones(endPoint)
    C6=2*np.ones(endPoint)
    C7=2*np.ones(endPoint)
    C8=2*np.ones(endPoint)
    
    print"\n============ 10 arrays (dim=%s) ================"%endPoint
    
    anaRes = anaResCalc(endPoint,2**8)
    print "anaRes: %s"%(anaRes)
    
    
    #numpy
    start = T()
    resNumpy = prodSumNumpy(A,B,C1,C2,C3,C4,C5,C6,C7,C8)
    end = T() - start
    print"numpy: %s (time: %1.2e)"%(resNumpy,end)
    
    #einsum
    start = T()
    resEinsum = prodSumEinsum(A,B,C1,C2,C3,C4,C5,C6,C7,C8)
    end = T() - start
    print"einsum: %s (time: %1.2e)"%(resEinsum,end)
