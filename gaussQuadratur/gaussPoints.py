"""
test env for gaussPoint calculation
"""
import numpy as np

class gaussPoints(object):
    def __init__(self,deg=5):
        self.__deg = deg
        self.resetPoints()
    
    def resetPoints(self):
        self.points, self.weights = np.polynomial.legendre.leggauss(self.__deg)
        self.bounds = (-1.0,1.0)
    
    def transform(self,low,up):
        if not(self.bounds == (-1.0,1.0)):
            self.resetPoints()
        self.points = (up-low)/2.0*self.points + (up+low)/2.0
        self.weights = (up-low)/2.0*self.weights
        self.bounds = (low,up)

if __name__=='__main__':
    degree = 5
    newBounds = (3.0,5.0)
    newnewBounds = (-10,10)
    
    gaussObj = gaussPoints(degree)
    print "bound: %s"%(str(gaussObj.bounds))
    print "points: %s"%(str(gaussObj.points))
    
    print "TRANSFORMATION 1"
    gaussObj.transform(*newBounds)
    print "bound: %s"%(str(gaussObj.bounds))
    print "points: %s"%(str(gaussObj.points))
    
    print "TRANSFORMATION 3"
    gaussObj.transform(*newnewBounds)
    print "bound: %s"%(str(gaussObj.bounds))
    print "points: %s"%(str(gaussObj.points))
    
    print "TRANSFORMATION 4 (retundance)"
    gaussObj.transform(*newnewBounds)
    print "bound: %s"%(str(gaussObj.bounds))
    print "points: %s"%(str(gaussObj.points))
