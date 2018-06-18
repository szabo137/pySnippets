"""
test env for some example values of internal integrals
"""
import unittest
import numpy as np
import sftrident.pulseLib as pulse



class mytest(unittest.TestCase):
    def setUp(self):
        self.__numInt = pulse.getPulse('cos')[1]
        self.__anaInt = pulse.getPulse('cos',['analytic'])[1]
        self.exampleDPhi = 5.0
        self.exampleArgs = np.array([-5.0,-2.5,0.0,2.5,5.0])
        self.exampleRes1 = np.array([np.pi**2*np.sin(5.0)/(50-2*np.pi**2),-0.715099,0.0,0.715099,-np.pi**2*np.sin(5.0)/(50-2*np.pi**2)])
        self.exampleRes2 = np.array([1.4187,1.41606,0.0,1.41606,1.4187])
        self.exampleRes3 = np.array([-0.933857,-0.813893,0.0,0.813893,0.933857])
        self.exampleRes4 = np.array([-0.941143,-0.919382,0.0,0.919382,0.941143])
        self.eps = 1e-4

    def test_AnaVals1(self):
        """
        tests the analytic values for some example points: integral1
        """
        self.anaRes1 = self.__anaInt(0)(self.exampleArgs,self.exampleDPhi)
        for i, item in enumerate(self.exampleArgs):
            self.assertLessEqual(abs(self.anaRes1[i] - self.exampleRes1[i]),self.eps)

    def test_AnaVals2(self):
        """
        tests the analytic values for some example points: integral2
        """
        self.anaRes2 = self.__anaInt(1)(self.exampleArgs,self.exampleDPhi)
        for i, item in enumerate(self.exampleArgs):
            self.assertLessEqual(abs(self.anaRes2[i] - self.exampleRes2[i]),self.eps)

    def test_AnaVals3(self):
        """
        tests the analytic values for some example points: integral3
        """
        self.anaRes3 = self.__anaInt(2)(self.exampleArgs,self.exampleDPhi)
        for i, item in enumerate(self.exampleArgs):
            self.assertLessEqual(abs(self.anaRes3[i] - self.exampleRes3[i]),self.eps)

    def test_AnaVals4(self):
        """
        tests the analytic values for some example points: integral4
        """
        self.anaRes4 = self.__anaInt(3)(self.exampleArgs,self.exampleDPhi)
        for i, item in enumerate(self.exampleArgs):
            self.res = self.anaRes4[i]
            self.assertLessEqual(abs(self.anaRes4[i] - self.exampleRes4[i]),self.eps)

    def test_numVals1(self):
        """
        tests the numeric values for some example points: integral1
        """
        self.numRes1 = self.__numInt(0)(self.exampleArgs,self.exampleDPhi)
        for i, item in enumerate(self.exampleArgs):
            self.assertLessEqual(abs(self.numRes1[i] - self.exampleRes1[i]),self.eps)

    def test_numVals2(self):
        """
        tests the numeric values for some example points: integral2
        """
        self.numRes2 = self.__numInt(1)(self.exampleArgs,self.exampleDPhi)
        for i, item in enumerate(self.exampleArgs):
            self.res = self.numRes2[i]
            self.assertLessEqual(abs(self.numRes2[i] - self.exampleRes2[i]),self.eps)

    def test_numVals3(self):
        """
        tests the numeric values for some example points: integral3
        """
        self.numRes3 = self.__numInt(2)(self.exampleArgs,self.exampleDPhi)
        for i, item in enumerate(self.exampleArgs):
            self.assertLessEqual(abs(self.numRes3[i] - self.exampleRes3[i]),self.eps)

    def test_numVals4(self):
        """
        tests the numeric values for some example points: integral4
        """
        self.numRes4 = self.__numInt(3)(self.exampleArgs,self.exampleDPhi)
        for i, item in enumerate(self.exampleArgs):
            self.assertLessEqual(abs(self.numRes4[i] - self.exampleRes4[i]),self.eps)
