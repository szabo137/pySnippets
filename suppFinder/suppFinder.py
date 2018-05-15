"""
module to find the approx support of a function
"""
import numpy as np
import matplotlib.pylab as plt



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



def env(x,dx):
    return (np.cos(np.pi*x/(2.0*dx)))**2*(x<= dx)*(x>=-dx)

class tester(object):
    def __init__(self):
        self.resetCallCount()
    
    def resetCallCount(self):
        self.callCount = 0
    
    def __call__(self,x,dx):
        self.callCount +=len(x)
        return (x-1.21559541)*(x-10.12022857)*np.cos(np.cos(x))*(np.sin(x)*env(x,dx))**2
        



testFunc = tester()


DD = 25
eps=1e-6

args = np.linspace(0.0,DD*1.1,500)
vals = testFunc(args,DD)
print "callCount: %s"%testFunc.callCount
plt.plot(args,vals)

testFunc.resetCallCount()



lagPoints = np.polynomial.laguerre.laggauss(15)[0]
print lagPoints
valsLag = testFunc(lagPoints,DD)
plt.plot(lagPoints,valsLag,'o')
print valsLag
zeroPointsCanidate=np.where(np.abs(valsLag)<eps)[0]
print zeroPointsCanidate

def checkSerial(arr):
    return arr[0]==arr[1]-1 and arr[1]==arr[2]-1

def getTriple(arr):
    print "blub: %s"%(np.arange(len(arr)-2))
    for n in np.arange(len(arr)-2):
        temp = arr[n:n+3]
        print"t: %s"%temp
        if checkSerial(temp):
            return temp
    print "Warning: there is no serial triple"
    return np.zeros(3)

zeroPoints = getTriple(zeroPointsCanidate)
print zeroPoints
#next step
newBound = lagPoints[zeroPoints-1]
print newBound
#zeroTest
gausPkte = gaussPoints(10)
print newBound[-2:]
gausPkte.transform(*newBound[-2:])
newPoints = gausPkte.points
print newPoints
newVal = testFunc(newPoints,DD)
gausZeros = np.where(np.abs(newVal)<eps)[0]
print gausZeros # return true if len(...)==len(gausPkte)
print len(gausZeros)==len(gausPkte.points)
plt.plot(newPoints,newVal,'x')

#boundtest

print "number of pings: %s"%testFunc.callCount
plt.show()
