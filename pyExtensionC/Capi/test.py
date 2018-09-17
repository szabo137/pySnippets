import cos_module as cm
from math import cos
import numpy as np

def runClock():
    import time as T
    a=range(50000)
    aNP = np.arange(len(a))
    start = T.time()
    for el in aNP:
        cos(cos(cos(cos(el))))
    end=T.time() - start
    start = T.time()
    for el in aNP:
        cm.cos_func(el)
    end2=T.time() -start
    start = T.time()
    np.cos(np.cos(np.cos(np.cos(aNP))))
    end3=T.time() -start
    print "math: %1.2e"%end
    print "c: %1.2e"%end2
    print "np: %1.2e"%end3


if __name__=='__main__':
    runClock()
