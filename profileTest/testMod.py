"""
module to test
"""


import numpy as np
def count():
    from math import sqrt
    for x in range(10**5):
        sqrt(x)

def count2():
    from math import sqrt
    for x in np.arange(10**5):
        sqrt(x)

def count3():
    x=np.arange(1e5)
    np.sqrt(x)

def runAll():
    count()
    count2()
    count3()


if __name__=='__main__':
    runAll()
