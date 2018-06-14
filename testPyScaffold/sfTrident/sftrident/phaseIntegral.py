"""
contains class to evaluate the phase integrals
"""
import numpy as np
#from sfTrident.phaseInt.pulseLib.cosPulse.pyEnvelope import envelope
#from sfTrident.phaseInt.pulseLib.cosPulse.analyticInt import analyticInternalInt as anaInt
import pulseLib
from scipy.special import airy

maxPhotoNum = 50.0 #config??


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



def photoAmp(item,envelope,xi=0.0,dPhi=10):
    """
    returns the photo amplitude as a function of phi
    """
    if item==3:
        return lambda x: (photoAmp(1,envelope,xi,dPhi)(x))**2 + (photoAmp(2,envelope,xi,dPhi)(x))**2
    else:
        return lambda x: np.cos(xi)*np.cos(x)*envelope(x,dPhi)*(item==1) + np.sin(xi)*np.sin(x)*envelope(x,dPhi)*(item==2)

def prodSum(*arrays):
    """
    calculation of the product of several 1-D arrays elementwise and summation of the products

    uses numpy
    """
    return np.sum(np.prod(arrays,axis=0))

class phaseIntegral(object):
    def __init__(self,config = {'xi':0.0,'dPhi':10,'envelope':'cos','pulseOpt':[]},deg=5):
        self.__xi = config['xi']
        self.__dPhi = config['dPhi']
        self.__deg = deg
        self.__gauss = gaussPoints(self.__deg)
        self.__gauss.transform(-self.__dPhi,self.__dPhi)

        self.__pulse = pulseLib.getPulse(config['envelope'],config['pulseOpt'])
        self.__envelope = self.__pulse[0]
        self.__internalInt = self.__pulse[1]

        self.setClean()
        self.evalPhiDependence()

    def setClean(self):
        """
        sets the used attributes to zero
        """
        self.__photoNum = 0.0
        self.__alphas = [0.0,0.0,0.0]
        self.B1=None #0.0 + 1j*0.0
        self.B2=None #0.0 + 1j*0.0
        self.B3=None #0.0 + 1j*0.0
        self.sTerm=np.zeros(self.__deg)

    def evalPhotoAmp(self,item,phiArray):
        """
        returns the evaluation of photoamp at the given gauss points
        """
        return photoAmp(item,self.__envelope,self.__xi,self.__dPhi)(phiArray)

    def evalBeta1(self,phiArray):
        """
        returns the evaluation of beta1 at the given gauss points
        """
        res = np.zeros(self.__deg)
        if np.cos(self.__xi)==0:
            return res
        else:
            return np.cos(self.__xi)*(self.__internalInt(0)(phiArray,self.__dPhi))

    def evalBeta2(self,phiArray):
        """
        returns the evaluation of beta2 at the given gauss points
        """
        res = np.zeros(self.__deg)
        if np.sin(self.__xi)==0:
            return res
        else:
            return np.sin(self.__xi)*(self.__internalInt(1)(phiArray,self.__dPhi))

    def evalBeta3(self,phiArray):
        """
        returns the evaluation of beta3 at the given gauss points
        """
        res = np.zeros(self.__deg)
        if np.sin(self.__xi)==0:
            return (np.cos(self.__xi))**2*(self.__internalInt(2)(phiArray,self.__dPhi))
        elif np.cos(self.__xi)==0:
            return (np.sin(self.__xi))**2*(self.__internalInt(3)(phiArray,self.__dPhi))
        else:
            return (np.cos(self.__xi))**2*(self.__internalInt(2)(phiArray,self.__dPhi)) + (np.sin(self.__xi))**2*(self.__internalInt(3)(phiArray,self.__dPhi))


    def evalPhiDependence(self):
        """
        evaluates all phi dependent quantities and set self.__array per eval
        """
        points = self.__gauss.points
        self.photoAmp1 = self.evalPhotoAmp(1,points)
        self.photoAmp2 = self.evalPhotoAmp(2,points)
        self.photoAmp3 = self.evalPhotoAmp(3,points)
        self.beta1 = self.evalBeta1(points)
        self.beta2 = self.evalBeta2(points)
        self.beta3 = self.evalBeta3(points)
        return 0

    def setPhotoNum(self,photoNum):
        self.__photoNum = photoNum
        #self.sTerm = np.exp(1j*self.__photoNum*self.__gauss.points)
        if np.abs(photoNum)>=maxPhotoNum:
            self.sTerm = np.zeros(self.__deg)
        else:
            self.sTerm = np.exp(1j*self.__photoNum*self.__gauss.points)
        return 0

    def setKin(self,alphas):
        self.__alphas = alphas
        self.A1 = np.exp(1j*self.__alphas[0]*self.beta1)
        self.A2 = np.exp(1j*self.__alphas[1]*self.beta2)
        self.A3 = np.exp(1j*self.__alphas[2]*self.beta3)
        return 0

    def getB1(self):
        self.B1 = prodSum(self.__gauss.weights,self.photoAmp1,self.sTerm,self.A1,self.A2,self.A3)
        #print self.B1
        return self.B1

    def getB2(self):
        self.B2 = prodSum(self.__gauss.weights,self.photoAmp2,self.sTerm,self.A1,self.A2,self.A3)
        return self.B2

    def getB3(self):
        self.B3 = prodSum(self.__gauss.weights,self.photoAmp3,self.sTerm,self.A1,self.A2,self.A3)
        return self.B3

    def getB0(self):
        if self.B1==None:
            self.getB1()

        if self.B2==None:
            self.getB2()

        if self.B3==None:
            self.getB3()

        #if self.__photoNum==0:
        #    raise ValueError("You need to set the photoNumber first:\n\n\t <phaseIntClass>.setPhotoNum(s!=0)")
        self.B0 = -1.0/self.__photoNum*prodSum(self.__alphas,[self.B1,self.B2,self.B3])
        return self.B0

    def B1func(self,photoNum):
        sTerm = np.where(photoNum<maxPhotoNum,np.exp(1j*photoNum*self.__gauss.points),0.0)
        return prodSum(self.__gauss.weights,self.photoAmp1,sTerm,self.A1,self.A2,self.A3)

    def B2func(self,photoNum):
        sTerm = np.where(photoNum<maxPhotoNum,np.exp(1j*photoNum*self.__gauss.points),0.0)
        return prodSum(self.__gauss.weights,self.photoAmp2,sTerm,self.A1,self.A2,self.A3)

    def B3func(self,photoNum):
        sTerm = np.where(photoNum<maxPhotoNum,np.exp(1j*photoNum*self.__gauss.points),0.0)
        return prodSum(self.__gauss.weights,self.photoAmp3,sTerm,self.A1,self.A2,self.A3)

    def B0func(self,photoNum):
        #self.setPhotoNum(photoNum)
        #self.getB1()
        #self.getB2()
        #self.getB3()
        B1 = self.B1func(photoNum)
        B2 = self.B2func(photoNum)
        B3 = self.B3func(photoNum)
        #return -1.0/photoNum*prodSum(self.__alphas,[self.B1,self.B2,self.B3])
        return -1.0/photoNum*prodSum(self.__alphas,[B1,B2,B3])


    def getGauss(self):
        """
        for testing purpose only

        returns the gauss object
        """
        return self.__gauss

class phaseIntegralCCF(object):
    def __init__(self,config = {},deg=0):
        pass

    def __b(self,e,c1,c2):
        return 2.0*np.pi*np.exp(1j*e)/((3.0*c2)**(1.0/3.0))

    def __eta(self,r,c1,c2):
        return -r*c1/3.0/c2 + 2.0*c1**3/(27*c2**2)

    def __mu(self,r,c1,c2):
        return (r-c1**2/(3.0*c2))/((3.0*c2)**(1.0/3.0))

    def setPhotoNum(self,photoNum):
        self.__photoNum = photoNum
        return 0

    def setKin(self,alphas):
        """
        only alpha1 and alpha3 are used in ccf, for technical propose: alpha input is [alpha1,alpha2,alpha3]

        todo:
            -   improve this!
        """
        self.__alphas = alphas
        self.__c1,self.__c2 = self.__alphas[0]/2.0,self.__alphas[2]/3.0
        return 0

    def B0func(self,r):
        et = self.__eta(r,self.__c1,self.__c2)
        prefac = self.__b(et,self.__c1,self.__c2)
        arg = self.__mu(r,self.__c1,self.__c2)
        return prefac*airy(arg)[0]

    def B1func(self,r):
        et = self.__eta(r,self.__c1,self.__c2)
        prefac = self.__b(et,self.__c1,self.__c2)
        arg = self.__mu(r,self.__c1,self.__c2)
        fac1 = self.__c1/3.0/self.__c2*airy(arg)[0]
        fac2 =1j/((3.0*self.__c2)**(1.0/3.0))*airy(arg)[1]
        return -prefac*(fac1 + fac2)

    def B3func(self,r):
        et = self.__eta(r,self.__c1,self.__c2)
        prefac = self.__b(et,self.__c1,self.__c2)
        arg = self.__mu(r,self.__c1,self.__c2)
        fac1 = ((self.__c1/3.0/self.__c2)**2 - arg/((3.0*self.__c2)**(2.0/3.0)))*airy(arg)[0]
        fac2 =1j*2.0*self.__c1/((3.0*self.__c2)**(4.0/3.0))*airy(arg)[1]
        return prefac*(fac1 + fac2)

    def B2func(self,r):
        """
        there is no B2 in ccf, but nessesay for following code
        todo:
            -   improve this!
        """
        if np.isscalar(r):
            temp = np.array(r)
        else:
            temp = r
        return np.zeros(r.shape)


    def getB1(self):
        self.B1 = self.B1func(self.__photoNum)
        return self.B1

    def getB2(self):
        self.B2 = self.B1func(self.__photoNum)
        return self.B2

    def getB3(self):
        self.B3 = self.B1func(self.__photoNum)
        return self.B3

    def getB0(self):
        self.B0 = self.B1func(self.__photoNum)
        return self.B0
