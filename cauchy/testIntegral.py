import numpy as np
from scipy import integrate, special
import time


def func(x):
    #print "ping"
    #time.sleep(0.01)
    return np.exp(100*x)


start = time.time()
resQUAD = integrate.quad(func, -1, 1, weight='cauchy', wvar=0)
end=time.time() - start

print"resQUAD: %s (time: %1.2e)"%(resQUAD,end)


# Check against known result
#print(2*special.sici(1)[0])

def prodSum(*arrays):
    """
    calculation of the product of several 1-D arrays elementwise and summation of the products

    uses numpy
    """
    return np.sum(np.prod(arrays,axis=0))


class gaussPoints(object):
    def __init__(self,deg=5):
        self.__deg = deg
        self.resetPoints()

    def resetPoints(self):
        self.points, self.weights = np.polynomial.legendre.leggauss(self.__deg)
        #self.points, self.weights = np.polynomial.chebyshev.chebgauss(self.__deg)
        self.bounds = (-1.0,1.0)

    def transform(self,low,up):
        if not(self.bounds == (-1.0,1.0)):
            self.resetPoints()
        self.points = (up-low)/2.0*self.points + (up+low)/2.0
        self.weights = (up-low)/2.0*self.weights
        self.bounds = (low,up)





start = time.time()
gaussObj=gaussPoints(100)
gaussObj.transform(0.0,1.0)

evalFunc=lambda x: (func(x)-func(-x))/x

vals = evalFunc(gaussObj.points)

resGauss= prodSum(gaussObj.weights,vals)
end=time.time() - start
print"resGauss: %s (time: %1.2e)"%(resGauss,end)
