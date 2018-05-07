import numpy as np
from scipy import integrate, special
import time
import math


def func(x):
    #print "ping"
    #time.sleep(0.01)
    return np.array([np.sin(x) for item in np.arange(16)])


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




def prodSum(*arrays):
    """
    calculation of the product of several 1-D arrays elementwise and summation of the products

    uses numpy
    """
    return np.sum(np.prod(arrays,axis=0))


class gaussPoints(object):
    def __init__(self,deg=5,mode="gaussLeg"):
        self.__deg = deg
        self.mode=mode
        self.__setPoints()
        self.__resetPoints()

    def __setPoints(self):
        if self.mode=="gaussLeg":
            polynomial=np.polynomial.legendre.leggauss
        elif self.mode=="gaussCheb":
            polynomial=np.polynomial.chebyshev.chebgauss
        elif self.mode=="gaussLag":
            polynomial=numpy.polynomial.laguerre.laggauss
        else:
            raise ValueError("<%s> is not a mode of gaussPoints!")

        self.__initPoints, self.__initWeights=polynomial(self.__deg)

    def __resetPoints(self):
        self.points, self.weights = self.__initPoints,self.__initWeights
        self.bounds = (-1.0,1.0)

    def transform(self,low,up):
        self.points = (up-low)/2.0*self.__initPoints + (up+low)/2.0
        self.weights = (up-low)/2.0*self.__initWeights
        self.bounds = (low,up)






gaussObj=gaussPoints(1000)
evalFunc=lambda x: (func(x)-func(-x))/x


start = time.time()
resGauss= integCauchy(func,1000.0,gaussObj)
end=time.time() - start
print"resGauss: %s (time: %1.2e)"%(resGauss,end)
print"rel Err: %s"%(np.abs(np.pi-resGauss)/np.pi)



start = time.time()
map(func,np.linspace(1,10,100))
end=time.time() - start
print"funcEval: %s rounds (time: %1.2e)"%(100,end)
