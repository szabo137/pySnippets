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

if __name__=='__main__':
    from time import time as T
    endPoint = 100000
    a=np.arange(endPoint)
    b=2.0*np.arange(endPoint)
    anaRes = (endPoint-1)*endPoint*(2*(endPoint-1) + 1)/3.0
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
