"""
tests how fast lambda vs partial vs exec is
"""
from timer import timer
import numpy as np
import time

def test(x):
    time.sleep(1e-4)
    return np.cos(x)


def pyCall(x,a,b,c):
    return x*a**3 + x**2*np.sqrt(b**2+1.0) + test(c)

@timer
def callWrapper(**kw):
    for xx in np.arange(1000):
        pyCall(xx,3,4,5)


func = lambda x: pyCall(x,3,4,5)

@timer
def lambdaWrapper(**kw):
    tempFunc = func
    for xx in np.arange(1000):
        tempFunc(xx)

class pyClass(object):
    def __init__(self,a,b,c):
        self.__args = [a,b,c]
        self.__preCalc()

    def __preCalc(self):
        self.__first = self.__args[0]**3
        self.__sec = np.sqrt(self.__args[1]**2+1.0)
        self.__third = test(self.__args[2])

    def eval(self,x):
        return self.__first*x + self.__sec*x**2 + self.__third


pyObj = pyClass(3,4,5)

@timer
def classWrapper():
    for xx in np.arange(1000):
        pyObj.eval(xx)



if __name__=='__main__':
    callWrapper()
    lambdaWrapper()
    classWrapper()
