"""
module to test own gaussian quadratur
"""

import numpy as np
from scipy.integrate import quad
import time

deg = 24
start = time.time()
gaussPoints, weights = np.polynomial.legendre.leggauss(deg)
end1 = time.time() - start


print"gauss points/weights: %s (time: %1.2e)"%(len(gaussPoints),end1)

def transformGauss(lowBound,upBound):
    return ((upBound-lowBound)/2.0*gaussPoints + (upBound+lowBound)/2.0,(upBound-lowBound)/2.0*weights)

start = time.time()
newGauss,newWeights = transformGauss(2,3)
end2 = time.time() - start

print"new gauss points/weights: %s (time:%1.2e)"%(len(newGauss),end2)

def func(x):
    return np.exp(1j*x)/x

def func1(x):
    return np.cos(x)/x

def func2(x):
    return np.sin(x)/x


start = time.time()
res = quad(func1,2,3)[0] + 1j*quad(func2,2,3)[0]
end3 = time.time() - start
print"quad: %s (time: %1.2e)"%(res,end3)

start = time.time()
funcPoints = func(newGauss)
end4 = time.time() - start

start = time.time()
integrand = newWeights*funcPoints
end5 = time.time() - start

start = time.time()
gaussRes = np.sum(integrand)
end6 = time.time() - start

print"mygauss: %s (full time: %1.2e, func: %1.2e, prod: %1.2e, sum: %1.2e)"%(gaussRes,end1+end2+end4+end5+end6,end4,end5,end6)


def buildGaussPoints(deg,lowBound=-1.0,upBound=1.0):
    gaussPoints,weights = np.polynomial.legendre.leggauss(deg)
    return ((upBound-lowBound)/2.0*gaussPoints + (upBound+lowBound)/2.0,(upBound-lowBound)/2.0*weights)

def gaussInteg(func,gaussPonts,weights):
    funcPoints = func(gaussPonts)
    integrand = weights*funcPoints
    return np.sum(integrand)

p,w = buildGaussPoints(5,2,3)
print gaussInteg(func,p,w)
