"""
here is the test
"""
import rootDir.mod.moduleToTest as mod
import numpy as np

import unittest

class mytest(unittest.TestCase):
    def setUp(self):
        self.func=mod.funcSQUARE
    def test1(self):
        self.assertEqual(self.func(2),4)
    

if __name__=='__main__':
    unittest.main()
