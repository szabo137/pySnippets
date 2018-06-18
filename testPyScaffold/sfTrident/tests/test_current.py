"""
test enviroment for current class

todo:
    -   implement test for J0 and J3 against tobias?!
"""
import unittest
import numpy as np
#from sfTrident.kinematics.kinClass import kinClass
#from sfTrident.kinematics.lab import getPhysArea
import sftrident.qft as qft

from sftrident.current import currentClass

import itertools

import tobiasLib.currents as tobias
m=1.0
a0=0.5
#xi=np.pi/4.0
theta_ph      = 0.
phi_ph        = 0.
#tobias sample momenta
pPos = qft.MinkowskiVector([ 1.41518852,  1. ,         0. ,         0.05252185])
pEl = qft.MinkowskiVector([ 1.41518852 ,-1.   ,       0.      ,    0.05252185])
pKlaser = qft.MinkowskiVector([ 1.46771037 , 0.   ,       0.    ,     -1.46771037])
pKphoton = qft.MinkowskiVector([ 1.46771037 , 0.    ,      0.    ,      1.46771037])

def eps_photon(l):
    """
    single photon polariation
    """
    return (l==1)*qft.MinkowskiVector([0,np.cos(theta_ph)*np.cos(phi_ph),np.cos(theta_ph)*np.sin(phi_ph),-np.sin(theta_ph)])\
           +(l==2)*qft.MinkowskiVector( [0,-np.sin(phi_ph),np.cos(phi_ph),0])

def eps_laser(l):
    """
    laser polarisation
    """
    return (l==1)*qft.MinkowskiVector( [0,1,0,0] ) + (l==2)*qft.MinkowskiVector( [0,0,1,0] )

def buildMom(E,cTh):
    if not(isinstance(E,np.ndarray) and isinstance(cTh,np.ndarray)):
        if np.isscalar(E):
            E=[E]
        if np.isscalar(cTh):
            cTh=[cTh]
        E=np.array(E)
        cTh = np.array(cTh)
    rho = np.sqrt(E**2 - np.ones(E.shape))
    return qft.MinkowskiVector([E,rho*np.sqrt(np.ones(cTh.shape)-cTh*cTh),0,cTh*rho])

def buildRandomMom():
    E = np.random.uniform(1.1,3.0,100)
    cTh = np.random.uniform(-1.0,1.0,100)
    momTemp = buildMom(E,cTh)
    return map(qft.MinkowskiVector,(momTemp.data).T)



class mytest(unittest.TestCase):
    def setUp(self):
        #self.kinObj=kinClass()
        #self.kinParaSet = getSetOfKinPara()
        self.xi =np.pi/4.0
        self.currentObj = currentClass({'mass':1.0,'a0':a0,'xi':self.xi})
        self.randomMom = buildRandomMom()
        self.eps=1e-12
        self.sampleMom = np.array([pPos,pEl,pKlaser,pKphoton])
        self.Vsample = np.array([qft.SpinorV ((self.sampleMom[0], m),s1) for s1 in [1,2]])
        self.UbarSample  = np.array([qft.SpinorUBar((self.sampleMom[1],m),s1) for s1 in [1,2]])
        self.pol1 = eps_laser(1)*np.cos(self.xi)-1j*eps_laser(2)*np.sin(self.xi) #for J1
        self.pol2 = eps_laser(1)*np.cos(self.xi)+1j*eps_laser(2)*np.sin(self.xi) #for J2
        self.polPhoto = np.array([eps_photon(l) for l in [1,2]])

    def test1(self):
        """
        J0 test 1: J0(ubarP,ubarP)*P == 2.0 for same spin, as well
        J0(ubarP,ubarP)*P == 0.0 for different spin
        """
        #print "No. of Momenta: %s"%(len(self.randomMom))
        for el in self.randomMom:
            U=np.array([qft.SpinorU((el,1.0),s) for s in [1,2]])
            Ubar=np.array([qft.SpinorUBar((el,1.0),s) for s in [1,2]])
            self.res1 = self.currentObj.generalJ0(Ubar,U)
            for el2 in self.res1[:2]:
                self.test = np.abs(el*el2-2.0)
                self.assertLessEqual(self.test,self.eps)
            for el2 in self.res1[2:]:
                self.test = np.abs(el*el2)
                self.assertLessEqual(self.test,self.eps)

    def testGeneralJ1pol1(self):
        """
        tests generalJ1 against tobias with several combinations of spin
        on one point and photoPol = 0
        """
        #tobias currents
        self.J1tobias00 = tobias.J1plus(self.UbarSample[0],self.Vsample[0],self.polPhoto[0],pKlaser,pEl,pPos)
        self.J1tobias11 = tobias.J1plus(self.UbarSample[1],self.Vsample[1],self.polPhoto[0],pKlaser,pEl,pPos)
        self.J1tobias01 = tobias.J1plus(self.UbarSample[0],self.Vsample[1],self.polPhoto[0],pKlaser,pEl,pPos)
        self.J1tobias10 = tobias.J1plus(self.UbarSample[1],self.Vsample[0],self.polPhoto[0],pKlaser,pEl,pPos)

        # my currents
        self.testJ1 = self.currentObj.generalJ12(self.sampleMom[2],self.sampleMom[1],self.sampleMom[0],self.UbarSample,self.Vsample,self.pol1)*self.polPhoto[0]

        #test
        self.assertLessEqual(np.abs(self.testJ1[0]/2.0-self.J1tobias00),self.eps)
        self.assertLessEqual(np.abs(self.testJ1[1]/2.0-self.J1tobias11),self.eps)
        self.assertLessEqual(np.abs(self.testJ1[2]/2.0-self.J1tobias01),self.eps)
        self.assertLessEqual(np.abs(self.testJ1[3]/2.0-self.J1tobias10),self.eps)

    def testGeneralJ1pol2(self):
        """
        tests generalJ1 against tobias with several combinations of spin
        on one point and photoPol = 1
        """
        #tobias currents
        self.J1tobias00 = tobias.J1plus(self.UbarSample[0],self.Vsample[0],self.polPhoto[1],pKlaser,pEl,pPos)
        self.J1tobias11 = tobias.J1plus(self.UbarSample[1],self.Vsample[1],self.polPhoto[1],pKlaser,pEl,pPos)
        self.J1tobias01 = tobias.J1plus(self.UbarSample[0],self.Vsample[1],self.polPhoto[1],pKlaser,pEl,pPos)
        self.J1tobias10 = tobias.J1plus(self.UbarSample[1],self.Vsample[0],self.polPhoto[1],pKlaser,pEl,pPos)
        self.__tobType = type(self.J1tobias00)
        # my currents
        self.testJ1 = self.currentObj.generalJ12(self.sampleMom[2],self.sampleMom[1],self.sampleMom[0],self.UbarSample,self.Vsample,self.pol1)*self.polPhoto[1]
        self.__mytype = type(self.testJ1)
        #test
        self.assertLessEqual(np.abs(self.testJ1[0]/2.0-self.J1tobias00),self.eps)
        self.assertLessEqual(np.abs(self.testJ1[1]/2.0-self.J1tobias11),self.eps)
        self.assertLessEqual(np.abs(self.testJ1[2]/2.0-self.J1tobias01),self.eps)
        self.assertLessEqual(np.abs(self.testJ1[3]/2.0-self.J1tobias10),self.eps)

    def testJ1pol1(self):
        """
        tests J1 against tobias with several combinations of spin
        on one point
        """
        #tobias currents
        self.J1tobias00 = tobias.J1plus(self.UbarSample[0],self.Vsample[0],self.polPhoto[0],pKlaser,pEl,pPos)
        self.J1tobias11 = tobias.J1plus(self.UbarSample[1],self.Vsample[1],self.polPhoto[0],pKlaser,pEl,pPos)
        self.J1tobias01 = tobias.J1plus(self.UbarSample[0],self.Vsample[1],self.polPhoto[0],pKlaser,pEl,pPos)
        self.J1tobias10 = tobias.J1plus(self.UbarSample[1],self.Vsample[0],self.polPhoto[0],pKlaser,pEl,pPos)
        self.__tobType = type(self.J1tobias00)
        # my currents
        self.currentObj.setKin(self.sampleMom[2],self.sampleMom[1],self.sampleMom[0],self.UbarSample,self.Vsample)
        self.testJ1 = self.currentObj.evalJ1()*self.polPhoto[0]
        self.testJ2 = self.currentObj.evalJ2()*self.polPhoto[0]
        self.testJplus = 0.5*(np.cos(self.xi)*self.testJ1 - 1j*np.sin(self.xi)*self.testJ2)
        self.__mytype = type(self.testJ1)
        #test
        self.assertLessEqual(np.abs(self.testJplus[0]-self.J1tobias00),self.eps)
        self.assertLessEqual(np.abs(self.testJplus[1]-self.J1tobias11),self.eps)
        self.assertLessEqual(np.abs(self.testJplus[2]-self.J1tobias01),self.eps)
        self.assertLessEqual(np.abs(self.testJplus[3]-self.J1tobias10),self.eps)


    def testJ1pol2(self):
        """
        tests J1 against tobias with several combinations of spin
        on one point
        """
        #tobias currents
        self.J1tobias00 = tobias.J1plus(self.UbarSample[0],self.Vsample[0],self.polPhoto[1],pKlaser,pEl,pPos)
        self.J1tobias11 = tobias.J1plus(self.UbarSample[1],self.Vsample[1],self.polPhoto[1],pKlaser,pEl,pPos)
        self.J1tobias01 = tobias.J1plus(self.UbarSample[0],self.Vsample[1],self.polPhoto[1],pKlaser,pEl,pPos)
        self.J1tobias10 = tobias.J1plus(self.UbarSample[1],self.Vsample[0],self.polPhoto[1],pKlaser,pEl,pPos)
        self.__tobType = type(self.J1tobias00)
        # my currents
        self.currentObj.setKin(self.sampleMom[2],self.sampleMom[1],self.sampleMom[0],self.UbarSample,self.Vsample)
        self.testJ1 = self.currentObj.evalJ1()*self.polPhoto[1]
        self.testJ2 = self.currentObj.evalJ2()*self.polPhoto[1]
        self.testJplus = 0.5*(np.cos(self.xi)*self.testJ1 - 1j*np.sin(self.xi)*self.testJ2)
        self.__mytype = type(self.testJ1)
        #test
        self.assertLessEqual(np.abs(self.testJplus[0]-self.J1tobias00),self.eps)
        self.assertLessEqual(np.abs(self.testJplus[1]-self.J1tobias11),self.eps)
        self.assertLessEqual(np.abs(self.testJplus[2]-self.J1tobias01),self.eps)
        self.assertLessEqual(np.abs(self.testJplus[3]-self.J1tobias10),self.eps)

    def testJ2pol1(self):
        """
        tests J2 against tobias with several combinations of spin
        on one point
        """

        #tobias currents
        self.J1tobias00 = tobias.J1minus(self.UbarSample[0],self.Vsample[0],self.polPhoto[0],pKlaser,pEl,pPos)
        self.J1tobias11 = tobias.J1minus(self.UbarSample[1],self.Vsample[1],self.polPhoto[0],pKlaser,pEl,pPos)
        self.J1tobias01 = tobias.J1minus(self.UbarSample[0],self.Vsample[1],self.polPhoto[0],pKlaser,pEl,pPos)
        self.J1tobias10 = tobias.J1minus(self.UbarSample[1],self.Vsample[0],self.polPhoto[0],pKlaser,pEl,pPos)
        self.__tobType = type(self.J1tobias00)
        # my currents
        self.currentObj.setKin(self.sampleMom[2],self.sampleMom[1],self.sampleMom[0],self.UbarSample,self.Vsample)
        self.testJ1 = self.currentObj.evalJ1()*self.polPhoto[0]
        self.testJ2 = self.currentObj.evalJ2()*self.polPhoto[0]
        self.testJminus = 0.5*(np.cos(self.xi)*self.testJ1 + 1j*np.sin(self.xi)*self.testJ2)
        self.__mytype = type(self.testJ2)
        #test
        self.assertLessEqual(np.abs(self.testJminus[0]-self.J1tobias00),self.eps)
        self.assertLessEqual(np.abs(self.testJminus[1]-self.J1tobias11),self.eps)
        self.assertLessEqual(np.abs(self.testJminus[2]-self.J1tobias01),self.eps)
        self.assertLessEqual(np.abs(self.testJminus[3]-self.J1tobias10),self.eps)


    def testJ2pol2(self):
        """
        tests J2 against tobias with several combinations of spin
        on one point
        """

        #tobias currents
        self.J1tobias00 = tobias.J1minus(self.UbarSample[0],self.Vsample[0],self.polPhoto[1],pKlaser,pEl,pPos)
        self.J1tobias11 = tobias.J1minus(self.UbarSample[1],self.Vsample[1],self.polPhoto[1],pKlaser,pEl,pPos)
        self.J1tobias01 = tobias.J1minus(self.UbarSample[0],self.Vsample[1],self.polPhoto[1],pKlaser,pEl,pPos)
        self.J1tobias10 = tobias.J1minus(self.UbarSample[1],self.Vsample[0],self.polPhoto[1],pKlaser,pEl,pPos)
        self.__tobType = type(self.J1tobias00)
        # my currents
        self.currentObj.setKin(self.sampleMom[2],self.sampleMom[1],self.sampleMom[0],self.UbarSample,self.Vsample)
        self.testJ1 = self.currentObj.evalJ1()*self.polPhoto[1]
        self.testJ2 = self.currentObj.evalJ2()*self.polPhoto[1]
        self.testJminus = 0.5*(np.cos(self.xi)*self.testJ1 + 1j*np.sin(self.xi)*self.testJ2)
        self.__mytype = type(self.testJ2)
        #test
        self.assertLessEqual(np.abs(self.testJminus[0]-self.J1tobias00),self.eps)
        self.assertLessEqual(np.abs(self.testJminus[1]-self.J1tobias11),self.eps)
        self.assertLessEqual(np.abs(self.testJminus[2]-self.J1tobias01),self.eps)
        self.assertLessEqual(np.abs(self.testJminus[3]-self.J1tobias10),self.eps)

    def testJ3pol1(self):
        """
        tests J3 against tobias with several combinations of spin
        on one point
        """

        #tobias currents
        self.J3tobias00 = tobias.J02(self.UbarSample[0],self.Vsample[0],self.polPhoto[0],pKlaser,pEl,pPos)
        self.J3tobias11 = tobias.J02(self.UbarSample[1],self.Vsample[1],self.polPhoto[0],pKlaser,pEl,pPos)
        self.J3tobias01 = tobias.J02(self.UbarSample[0],self.Vsample[1],self.polPhoto[0],pKlaser,pEl,pPos)
        self.J3tobias10 = tobias.J02(self.UbarSample[1],self.Vsample[0],self.polPhoto[0],pKlaser,pEl,pPos)
        self.__tobType = type(self.J3tobias00)
        # my currents
        self.currentObj.setKin(self.sampleMom[2],self.sampleMom[1],self.sampleMom[0],self.UbarSample,self.Vsample)
        self.testJ3 = self.currentObj.evalJ3()*self.polPhoto[0]
        self.testJ02 = 0.5*self.testJ3
        self.__mytype = type(self.testJ3)
        #test
        self.assertLessEqual(np.abs(self.testJ02[0]-self.J3tobias00),self.eps)
        self.assertLessEqual(np.abs(self.testJ02[1]-self.J3tobias11),self.eps)
        self.assertLessEqual(np.abs(self.testJ02[2]-self.J3tobias01),self.eps)
        self.assertLessEqual(np.abs(self.testJ02[3]-self.J3tobias10),self.eps)

    def testJ3pol2(self):
        """
        tests J3 against tobias with several combinations of spin
        on one point
        """

        #tobias currents
        self.J3tobias00 = tobias.J02(self.UbarSample[0],self.Vsample[0],self.polPhoto[1],pKlaser,pEl,pPos)
        self.J3tobias11 = tobias.J02(self.UbarSample[1],self.Vsample[1],self.polPhoto[1],pKlaser,pEl,pPos)
        self.J3tobias01 = tobias.J02(self.UbarSample[0],self.Vsample[1],self.polPhoto[1],pKlaser,pEl,pPos)
        self.J3tobias10 = tobias.J02(self.UbarSample[1],self.Vsample[0],self.polPhoto[1],pKlaser,pEl,pPos)
        self.__tobType = type(self.J3tobias00)
        # my currents
        self.currentObj.setKin(self.sampleMom[2],self.sampleMom[1],self.sampleMom[0],self.UbarSample,self.Vsample)
        self.testJ3 = self.currentObj.evalJ3()*self.polPhoto[1]
        self.testJ02 = 0.5*self.testJ3
        self.__mytype = type(self.testJ3)
        #test
        self.assertLessEqual(np.abs(self.testJ02[0]-self.J3tobias00),self.eps)
        self.assertLessEqual(np.abs(self.testJ02[1]-self.J3tobias11),self.eps)
        self.assertLessEqual(np.abs(self.testJ02[2]-self.J3tobias01),self.eps)
        self.assertLessEqual(np.abs(self.testJ02[3]-self.J3tobias10),self.eps)
if __name__=='__main__':
    unittest.main()
