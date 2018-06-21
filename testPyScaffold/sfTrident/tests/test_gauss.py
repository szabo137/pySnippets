"""
test env for gaussIntegration
"""
import unittest
import numpy as np
import sftrident as sf

class mytest(unittest.TestCase):
    def setUp(self):
        self.degs = [2,5,10,100]
        self.newBounds=np.array([[3.0,5.0],[-10,10],[-37.5,-15.8],[-137,200]])
        self.examplePoints = np.array([-0.7745966692414834,0.0,0.7745966692414834])
        self.exampleWeights = np.array([0.555555555555555555555,0.88888888888888888888888,0.5555555555555555555555555])
        self.eps = 1e-14

    def test_output(self):
        """
        tests the length of output after initialisation
        """
        for el in self.degs:
            self.gauss = sf.math.gaussPoints(el)
            self.assertEqual(len(self.gauss.points),el)
            self.assertEqual(len(self.gauss.weights),el)
            self.assertEqual(self.gauss.bounds,(-1,1))

    def test_transfrom(self):
        """
        tests outputformat of transformed points and weights
        """
        for el in self.degs:
            self.gauss = sf.math.gaussPoints(el)
            for bounds in self.newBounds:
                self.testBound = bounds
                self.gauss.transform(*bounds)
                self.assertEqual(len(self.gauss.points),el)
                self.assertEqual(len(self.gauss.weights),el)
                self.assertEqual(self.gauss.bounds,tuple(bounds))

                #tests if all points in bounds
                for point in self.gauss.points:
                    self.assertGreaterEqual(point,self.gauss.bounds[0])
                    self.assertLessEqual(point,self.gauss.bounds[1])

    def test_retundance(self):
        """
        tests retundance of transformed points and weights
        """
        for el in self.degs:
            self.gauss = sf.math.gaussPoints(el)
            for bounds in self.newBounds:
                self.testBound = bounds
                self.gauss.transform(*bounds)
                self.gauss.transform(*bounds)
                self.assertEqual(self.gauss.bounds,tuple(bounds))
                self.gauss.transform(-1,1)
                self.assertEqual(self.gauss.bounds,(-1,1))

    def test_example(self):
        """
        tests values of degree = 3
        """
        self.gauss = sf.math.gaussPoints(3)
        for item,el in enumerate(self.gauss.points):
            self.assertLessEqual(abs(el - self.examplePoints[item]),self.eps)
        for item,el in enumerate(self.gauss.weights):
            self.assertLessEqual(abs(el - self.exampleWeights[item]),self.eps)
