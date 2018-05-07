import numpy as np
from scipy import integrate, special
import time
import math


def func(x):
    #print "ping"
    #time.sleep(0.01)
    return np.array([np.sin(x)*np.exp(1j*(x+np.sin(x)**2)) for item in np.arange(16)])


def integCauchy(fkt,bound,gauss=None,mode='quad'):
    """
    integrate fkt(x)/x over symetric interval [-bound,bound]
    """
    #bounds=gauss.bounds
    gauss.transform(0.0,bound)
    #test=np.ones(5)
    #test2=np.arange(5)
    tempFunc=lambda x: (fkt(x)-fkt(-x))/x
    #tempVals =(fkt(gauss.points) - fkt(-gauss.points))/gauss.points
    tempVals=np.array(map(tempFunc,gauss.points)).T
    #tempVals=np.map(tempFunc,gauss.points).T
    #tempVals =(fkt(test) - fkt(-test))/test
    #print tempVals
    #print np.array(tempVals2).T
    outTemp=(gauss.weights*tempVals)
    #outTemp=(test2*tempVals)
    #print outTemp
    #gauss.transform(*bounds)
    return np.sum(outTemp,axis=1)




startQUAD = time.time()
resQUAD = [integrate.quad(lambda x: func(x)[index], -1, 1, weight='cauchy', wvar=0,epsrel=1e-4) for index in np.arange(16)]
endQUAD=time.time() - startQUAD

print"resQUAD: %s (time: %1.2e)"%(resQUAD,endQUAD)


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






gaussObj=gaussPoints(200)
evalFunc=lambda x: (func(x)-func(-x))/x


startG = time.time()
resGauss= integCauchy(func,1.0,gaussObj)
endG=time.time() - startG
print"resGauss: %s (time: %1.2e)"%(resGauss,endG)

start = time.time()
map(func,np.linspace(1,10,100))
end=time.time() - start
print"funcEval: %s rounds (time: %1.2e)"%(100,end)
