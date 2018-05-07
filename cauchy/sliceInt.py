"""
integrate cauchy principal value with slicing the domain
"""
import numpy as np


def slicer(arr,x0=0,width=1.0):
    xMin = x0-width
    xMax = x0+width
    if (xMin in arr) and (xMax in arr):
        raise ValueError("The width is to wide! [%s,%s] is not a subset of [%s,%s]"%(xMin,xMax,arr[0],arr[-1]))
    return (arr[arr<(xMin)], arr[np.logical_and(arr>xMin,arr<xMax)],arr[arr>xMax])


if __name__=='__main__':
    a=np.linspace(-100,100,45)
    print"lower: %s\nmid: %s\nupper: %s"%slicer(a,10.0,15.4)
