"""
tests the behavior of parrays
"""

import numpy as np
import qft
import time

omega = qft.parray(np.arange(4))
omega2 = qft.parray(np.arange(3))

momArr = np.array([omega,0,0,omega])
momMks1 = qft.MinkowskiVector(momArr)
momArr = np.array([omega2,0,0,-omega2])
momMks2 = qft.MinkowskiVector(momArr)
momMks=momMks1*momMks2

facArr = np.array([1,2])
facParr = qft.parray(facArr)

try:
    start = time.time()
    res1=momArr*facArr
    end=time.time() - start
    print "momArr*facArr: %s %s"%(res1,end)
except:
    print "not possible!"

try:
    start = time.time()
    res2=momMks*facArr
    end=time.time() - start
    print "momMks*facArr: %s %s"%(res2,end)
except:
    print "not possible!"

try:
    start = time.time()
    res3=momArr*facParr
    end=time.time() - start
    print "momArr*facParr: %s %s"%(res3,end)
except:
    print "not possible!"

try:
    start = time.time()
    res4=np.outer(momMks,facParr)
    end=time.time() - start
    print "momMks*facParr: %s %s"%(res4,end)
except TypeError as err:
    print "not possible! %s"%err
