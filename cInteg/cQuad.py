"""
module to test ctypes and lowLevelCallable functions
"""

import os, ctypes
from scipy import integrate, LowLevelCallable

lib = ctypes.CDLL(os.path.abspath('testlib.so'))
lib.f.restype = ctypes.c_double
lib.f.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_void_p)

c=ctypes.c_double(1.0)
userData = ctypes.cast(ctypes.pointer(c), ctypes.c_void_p)
func = LowLevelCallable(lib.f, userData)

import time
t1=time.time()
res1 = integrate.nquad(func,[[0,10],[-10,0],[-1,1]])
t2=time.time()
print 'ctypes: %s (%s)'%(res1,t2-t1)

import numpy as np

def pyfunc(x,y,z):
    return 1.0 + np.exp(x)*np.sin(3.0*x) - y*z

t3=time.time()
res2 = integrate.nquad(pyfunc,[[0,10],[-10,0],[-1,1]])
t4=time.time()
print 'py: %s (%s)'%(res2,t4-t3)
