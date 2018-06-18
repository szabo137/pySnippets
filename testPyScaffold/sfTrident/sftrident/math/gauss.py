"""
The gauss module
================
This module provides gauss points and weights of several types.

Types
-----
 - gauss legendre (default)
 - gauss chebychev
 - gauss laguerre
"""
import numpy as np


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
