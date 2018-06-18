"""
test env to check if the analytic internal integrals and numeric evaluation give the same
"""
import unittest
import numpy as np
import sftrident.pulseLib as pulse



class mytest(unittest.TestCase):
    def setUp(self):
        self.__numInt = pulse.getPulse('cos')[1]
        self.__anaInt = pulse.getPulse('cos',['analytic'])[1]
        self.exampleDPhis = np.linspace(10,50,5)
        self.eps = 1e-13


    def test_Int1(self):
        """
        tests the integral No. 1:

        integrand1 = cos^2(phi'/(2 dPhi))*cos(phi')
        """
        for el in self.exampleDPhis:
            self.phiArray = np.linspace(-el,el,10)
            self.resNum = self.__numInt(0)(self.phiArray,el)
            self.resAna = self.__anaInt(0)(self.phiArray,el)
            for i,item in enumerate(self.resNum):
                self.assertLessEqual(abs(item - self.resAna[i]),self.eps)

    def test_Int2(self):
        """
        tests the integral No. 1:

        integrand1 = cos^2(phi'/(2 dPhi))*cos(phi')
        """
        for el in self.exampleDPhis:
            self.phiArray = np.linspace(-el,el,10)
            self.resNum = self.__numInt(1)(self.phiArray,el)
            self.resAna = self.__anaInt(1)(self.phiArray,el)
            for i,item in enumerate(self.resNum):
                self.assertLessEqual(abs(item - self.resAna[i]),self.eps)

    def test_Int3(self):
        """
        tests the integral No. 1:

        integrand1 = cos^2(phi'/(2 dPhi))*cos(phi')
        """
        for el in self.exampleDPhis:
            self.phiArray = np.linspace(-el,el,10)
            self.resNum = self.__numInt(2)(self.phiArray,el)
            self.resAna = self.__anaInt(2)(self.phiArray,el)
            for i,item in enumerate(self.resNum):
                self.assertLessEqual(abs(item - self.resAna[i]),self.eps)

    def test_Int4(self):
        """
        tests the integral No. 1:

        integrand1 = cos^2(phi'/(2 dPhi))*cos(phi')
        """
        for el in self.exampleDPhis:
            self.phiArray = np.linspace(-el,el,10)
            self.resNum = self.__numInt(3)(self.phiArray,el)
            self.resAna = self.__anaInt(3)(self.phiArray,el)
            for i,item in enumerate(self.resNum):
                self.assertLessEqual(abs(item - self.resAna[i]),self.eps)
