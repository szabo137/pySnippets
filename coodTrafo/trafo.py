"""
class to store several transformations -> to lightcone
"""
import numpy as np


def spherical2lightcone(E,cTh,phi,m=1.0):
    rho = np.sqrt(E**2 - m**2)*(E>=1)
    sTh = np.sqrt(1-cTh**2)
    return ((E-rho*cTh)/2.0,rho*sTh*np.cos(phi),rho*sTh*np.sin(phi))

def transverse2lightcone(y,pT,phi,m=1.0):
    return (0.5*np.exp(-y)*np.sqrt(pT**2 + m**2),pT*np.cos(phi),pT*np.sin(phi))


def jacTransverse():
    pass

def jacSpherical():
    pass

class toLightcone(object):
    def __init__(self,sourceCoord):
        self.source = sourceCoord
        if self.source=='spherical':
            self.__trafo = spherical2lightcone
            self.__jac = jacSpherical
        elif: self.source=='trans':
            self.__trafo = transverse2lightcone
            self.__jac = jacTransverse
        else:
            raise ValueError("The source coordinates <%s> are not known.")

    def getTrafo(self):
        return self.__trafo


    def getJac(self):
        return self.__jac
