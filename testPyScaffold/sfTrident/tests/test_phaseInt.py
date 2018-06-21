"""
test env for phaseIntClass
"""
import unittest
import numpy as np
import time

import sftrident.qft as qft
import sftrident as phase
from sftrident.kinematics import alpha

from tobiasLib.phaseInt import phaseIntFromTob as tobPI




class mytest(unittest.TestCase):
    def setUp(self):
        self.__tobB1,self.__tobB2,self.__tobB3,self.__tobB0,self.__sampleMomenta,para =tobPI()

        self.pPos = qft.MinkowskiVector(self.__sampleMomenta[0].data)
        #print "pPos: %s \n\t(type: %s)"%(self.pPos,type(self.pPos))
        #print qft.SpinorV((self.pPos,m),1)
        self.pEl = qft.MinkowskiVector(self.__sampleMomenta[1].data)
        #print "pEl: %s \n\t(type: %s)"%(self.pEl,type(self.pEl))
        self.pLaser = qft.MinkowskiVector(self.__sampleMomenta[2].data)
        #print "pLaser: %s \n\t(type: %s)"%(self.pLaser,type(self.pLaser))
        self.pPhoton = qft.MinkowskiVector(self.__sampleMomenta[3].data)
        #print "pPhoton: %s \n\t(type: %s)"%(self.pPhoton,type(self.pPhoton))

        self.bwMomenta = np.array([self.pLaser,self.pEl,self.pPos])


        self.__a0 = para[0]
        self.__xi = para[1]
        self.__dphi = para[2]
        self.__config = {'a0':self.__a0,'mass':1.0,'xi':self.__xi,'dPhi':self.__dphi,'envelope':'cos','pulseOpt':['analytic']}

        self.__photoNum =(self.pPos.minus() + self.pEl.minus() - self.pPhoton.minus())/self.pLaser.minus()
        self.alphas = np.array([alpha(self.bwMomenta,self.__config)(item) for item in [1,2,3]])
        #print self.__alphas

        self.__deg = 1000

        self.__eps = 1e-5#*(1.0 + 1j)

    def test_evalB1(self):
        """
        test evaluation of B_1 against tobias integrals
        """
        self.__phaseObj = phase.phaseIntegral(self.__config,deg=self.__deg)
        self.__phaseObj.setPhotoNum(self.__photoNum)
        self.__phaseObj.setKin(self.alphas)
        self.__res = self.__phaseObj.getB1()
        print"B1: %s"%self.__res
        print"tobB1: %s"%self.__tobB1
        self.assertLessEqual(np.abs(np.real(self.__res-self.__tobB1)),self.__eps)
        self.assertLessEqual(np.abs(np.imag(self.__res-self.__tobB1)),self.__eps)

    def test_evalB2(self):
        """
        test evaluation of B_2 against tobias integrals
        """

        self.__phaseObj = phase.phaseIntegral(self.__config,deg=self.__deg)
        self.__phaseObj.setPhotoNum(self.__photoNum)
        self.__phaseObj.setKin(self.alphas)
        self.__res = self.__phaseObj.getB2()
        print"B2: %s"%self.__res
        print"tobB2: %s"%self.__tobB2
        self.assertLessEqual(np.abs(np.real(self.__res-self.__tobB2)),self.__eps)
        self.assertLessEqual(np.abs(np.imag(self.__res-self.__tobB2)),self.__eps)


    def test_evalB3(self):
        """
        test evaluation of B_3 against tobias integrals
        """

        self.__phaseObj = phase.phaseIntegral(self.__config,deg=self.__deg)
        self.__phaseObj.setPhotoNum(self.__photoNum)
        self.__phaseObj.setKin(self.alphas)
        self.__res = self.__phaseObj.getB3()
        print"B3: %s"%self.__res
        print"tobB3: %s"%self.__tobB3
        self.assertLessEqual(np.abs(np.real(self.__res-self.__tobB3)),self.__eps)
        self.assertLessEqual(np.abs(np.imag(self.__res-self.__tobB3)),self.__eps)


    def test_evalB0(self):
        """
        test evaluation of B_0 against tobias integrals
        """

        self.__phaseObj = phase.phaseIntegral(self.__config,deg=self.__deg)
        self.__phaseObj.setPhotoNum(self.__photoNum)
        self.__phaseObj.setKin(self.alphas)
        self.__res = self.__phaseObj.getB0()
        print"B0: %s"%self.__res
        print"tobB0: %s"%self.__tobB0
        self.assertLessEqual(np.abs(np.real(self.__res-self.__tobB0)),self.__eps)
        self.assertLessEqual(np.abs(np.imag(self.__res-self.__tobB0)),self.__eps)

    def test_funcB1(self):
        """
        test methode funcB1 against evalB1
        """
        #eval B1
        self.__phaseObj = phase.phaseIntegral(self.__config,deg=self.__deg)
        self.__phaseObj.setKin(self.alphas)
        #print self.__a0
        for photoNum in np.linspace(-20,20,30):
            self.__phaseObj.setPhotoNum(photoNum)
            self.__resEval = self.__phaseObj.getB1()

            self.__resFunc = self.__phaseObj.B1func(photoNum)
            #print "eval: %s"%(self.__resEval)
            #print "func: %s"%(self.__resFunc)
            self.assertLessEqual(np.abs(np.real(self.__resEval-self.__resFunc)),self.__eps)
            self.assertLessEqual(np.abs(np.imag(self.__resEval-self.__resFunc)),self.__eps)

    def test_funcB2(self):
        """
        test methode funcB2 against evalB2
        """
        #eval B2
        self.__phaseObj = phase.phaseIntegral(self.__config,deg=self.__deg)
        self.__phaseObj.setKin(self.alphas)

        for photoNum in np.linspace(-20,20,30):
            start = time.time()
            self.__phaseObj.setPhotoNum(photoNum)
            self.__resEval = self.__phaseObj.getB2()
            end=time.time() - start
            #print "eval: %s (time: %1.2e)"%(self.__resEval,end)

            start = time.time()
            self.__resFunc = self.__phaseObj.B2func(photoNum)
            end=time.time() - start
            #print "func: %s (time: %1.2e)"%(self.__resFunc,end)
            self.assertLessEqual(np.abs(np.real(self.__resEval-self.__resFunc)),self.__eps)
            self.assertLessEqual(np.abs(np.imag(self.__resEval-self.__resFunc)),self.__eps)


    def test_funcB3(self):
        """
        test methode funcB3 against evalB3
        """
        #eval B2
        self.__phaseObj = phase.phaseIntegral(self.__config,deg=self.__deg)
        self.__phaseObj.setKin(self.alphas)

        for photoNum in np.linspace(-20,20,30):
            start = time.time()
            self.__phaseObj.setPhotoNum(photoNum)
            self.__resEval = self.__phaseObj.getB3()
            end=time.time() - start
            print "eval: %s (time: %1.2e)"%(self.__resEval,end)

            start = time.time()
            self.__resFunc = self.__phaseObj.B3func(photoNum)
            end=time.time() - start
            print "func: %s (time: %1.2e)"%(self.__resFunc,end)
            self.assertLessEqual(np.abs(np.real(self.__resEval-self.__resFunc)),self.__eps)
            self.assertLessEqual(np.abs(np.imag(self.__resEval-self.__resFunc)),self.__eps)
