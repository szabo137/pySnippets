"""
contains the integrands of the internal integration (ctypes callable)

todo:
 - low level callable work like have the pi in the envelope
"""
import numpy as np
from scipy.integrate import quad
import os
here = os.path.dirname(__file__)
os.system("gcc -shared -fPIC -o %s/cIntegrands.so %s/cIntegrands.c"%(here,here))
#print "compile done"

#from sfTrident.settings import laserConstants
from math import pi
#divided by pi because its missing in the cfunctions
#dPhi = laserConstants['dphi']/pi

from scipy import LowLevelCallable
import os, ctypes

lib = ctypes.CDLL(os.path.abspath('%s/cIntegrands.so'%here))
#lib = ctypes.CDLL(os.path.abspath('sfTrident/phaseInt/pulseLib/cosPulse/cIntegrands.so'))


def __I1(phi,dPhi):
    dPhi=dPhi/pi
    lib.integF1.restype = ctypes.c_double
    lib.integF1.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_void_p)
    c=ctypes.c_double(dPhi)
    userData = ctypes.cast(ctypes.pointer(c), ctypes.c_void_p)
    func1 = LowLevelCallable(lib.integF1, userData)
    if np.isscalar(phi):
        res = quad(func1,0,phi)[0]
    else:
        res = np.array([quad(func1,0,x)[0] for x in phi])
    return res

def __I2(phi,dPhi):
    dPhi=dPhi/pi
    lib.integF2.restype = ctypes.c_double
    lib.integF2.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_void_p)
    c=ctypes.c_double(dPhi)
    userData = ctypes.cast(ctypes.pointer(c), ctypes.c_void_p)
    func1 = LowLevelCallable(lib.integF2, userData)
    if np.isscalar(phi):
        res = quad(func1,0,phi)[0]
    else:
        res = np.array([quad(func1,0,x)[0] for x in phi])
    return res

def __I3(phi,dPhi):
    dPhi=dPhi/pi
    lib.integF3.restype = ctypes.c_double
    lib.integF3.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_void_p)
    c=ctypes.c_double(dPhi)
    userData = ctypes.cast(ctypes.pointer(c), ctypes.c_void_p)
    func1 = LowLevelCallable(lib.integF3, userData)
    if np.isscalar(phi):
        res = quad(func1,0,phi)[0]
    else:
        res = np.array([quad(func1,0,x)[0] for x in phi])
    return res

def __I4(phi,dPhi):
    dPhi=dPhi/pi
    lib.integF4.restype = ctypes.c_double
    lib.integF4.argtypes = (ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_void_p)
    c=ctypes.c_double(dPhi)
    userData = ctypes.cast(ctypes.pointer(c), ctypes.c_void_p)
    func1 = LowLevelCallable(lib.integF4, userData)
    if np.isscalar(phi):
        res = quad(func1,0,phi)[0]
    else:
        res = np.array([quad(func1,0,x)[0] for x in phi])
    return res

def getIntegrals(item):
    integrals = [__I1,__I2,__I3,__I4]
    return integrals[item]
