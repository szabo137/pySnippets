"""
module to test conservation law for momenta

todo:
    -   eval kinObj only once
    -   develop value test
"""
import unittest
import numpy as np
from sftrident.kinematics import kinClass
from sftrident.kinematics import modeDict
import sftrident.qft as qft

import itertools


def getSetOfKinPara(eps=1e-10):
    ss=np.array([4.0])
    p2x=np.linspace(-10,10,3)
    p2y=np.linspace(-10,10,3)
    p2minus =np.linspace(-0.1,0.5,3)
    p3x=np.linspace(-10,10,3)
    p3y=np.linspace(-10,10,3)
    p3minus = np.linspace(-0.1,0.5,3)
    kinPara=[t for t in itertools.product(ss,p2x,p2y,p2minus,p3x,p3y,p3minus)]
    return np.array(kinPara)


def calcConserve(momenta):
    """
    beware of the minus sign -p !!!!
    """
    #return np.sum(momenta[:2])-np.sum(momenta[3:])
    pt = momenta[2] + momenta[3] + momenta[4]
    rPluS = -(pt*pt - 1.0)/(2.0*(momenta[0]*momenta[1]))
    return pt + momenta[1] - rPluS*momenta[0]

def calcConserve2(momenta,r,s):
    return momenta[0]*(r+s) + momenta[1] - momenta[2] - momenta[3] - momenta[4]

class mytest(unittest.TestCase):
    def setUp(self):
        self.kinObj = kinClass(config={'a0':1.0,'mass':1.0})
        self.kinParaSet = getSetOfKinPara()
        self.conserve = calcConserve
        self.eps=1e-8
        self.momType=type(qft.MinkowskiVector([0,0,0,0]))
        self.spinorUType=type(qft.SpinorU(([1,0,0,1],0),1))
        self.spinorUBarType=type(qft.SpinorUBar(([1,0,0,1],0),1))
        self.spinorVType=type(qft.SpinorV(([1,0,0,1],0),1))
        self.spinorVBarType=type(qft.SpinorVBar(([1,0,0,1],0),1))

    def test_typetest(self):
        """
        typetest
        """
        for el in self.kinParaSet:
            self.kinObj.evalKin(el)
            mom=self.kinObj.getMomenta()
            spinors = self.kinObj.getSpinors()
            #print self.kinObj.physArea
            if self.kinObj.physArea:
                #mom type test
                for item in np.arange(len(mom)):
                    self.assertEqual(type(mom[item]),self.momType)
                #spinor type test
                self.assertEqual(type(spinors[0][0]),self.spinorUType)
                self.assertEqual(type(spinors[0][1]),self.spinorUType)
                for el3 in spinors[[1,2]]:
                    self.assertEqual(type(el3[0]),self.spinorUBarType)
                    self.assertEqual(type(el3[1]),self.spinorUBarType)
                self.assertEqual(type(spinors[3][0]),self.spinorVType)
                self.assertEqual(type(spinors[3][1]),self.spinorVType)

    def test_conservation(self):
        """
        tests the conversation of output momenta (direct)
        """
        for el in self.kinParaSet:
            self.kinObj.evalKin(el)
            #self.assertEqual(self.kinObj.physArea,True)
            if self.kinObj.physArea:
                mom=self.kinObj.getMomenta()
                self.assertEqual(len(mom),5)
                self.assertLessEqual(self.conserve(mom)()[0], self.eps)
                self.assertLessEqual(self.conserve(mom)()[1], self.eps)
                self.assertLessEqual(self.conserve(mom)()[2], self.eps)
                self.assertLessEqual(self.conserve(mom)()[3], self.eps)

    def test_modeDict(self):
        """
        check if selection in modeDict is correct
        """
        self.modeDict = modeDict
        self.momNames = np.array(['k','p','p1','p2','p3'])
        self.spinorNames = np.array(['u','ubar1','ubar2','v3'])
        self.alphaNames = np.array(['alphaBW','alphaC','alphaBWx','alphaCx'])
        #momTest
        for item in [0,1,2]:
            self.assertEqual(self.momNames[modeDict['bw'][0]][item],np.array(['k','p2','p3'])[item])
            self.assertEqual(self.momNames[modeDict['c'][0]][item],np.array(['k','p1','p'])[item])
            self.assertEqual(self.momNames[modeDict['bwx'][0]][item],np.array(['k','p1','p3'])[item])
            self.assertEqual(self.momNames[modeDict['cx'][0]][item],np.array(['k','p2','p'])[item])


if __name__=='__main__':
    unittest.main()
