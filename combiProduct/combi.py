"""
product of all combinations of entries of an array

the result is an len(a)xlen(b)-array
"""
import numpy as np
a=np.arange(100)
b=3*np.arange(100)

print a
print b

def combiProdFor(x,y):
    res = np.zeros((len(x),len(y)))
    for item,el in enumerate(x):
        for item2,el2 in enumerate(y):
            res[item][item2]=el*el2
    return res.reshape(len(x)*len(y),)

def combiProdNumpy1(x,y):
    return np.array([x[0]*y,x[1]*y,x[2]*y,x[3]*y]).reshape(len(x)*len(y),)

def combiProdNumpy2(x,y):
    return np.array([x[item]*y for item in np.arange(len(x))]).reshape(len(x)*len(y),)

def combiProdNumpy3(x,y):
    return np.outer(x,y).flatten()

import time as T
print "-"*30
start = T.time()
resFor = combiProdFor(a,b)
end = T.time() - start
print "for: %s \n time: %1.2e"%(resFor,end)

#print "-"*30
#start = T.time()
#resNumpy1 = combiProdNumpy1(a,b)
#end = T.time() - start
#print "Numpy1: %s \n time: %1.2e"%(resNumpy1,end)

print "-"*30
start = T.time()
resNumpy2 = combiProdNumpy2(a,b)
end = T.time() - start
print "Numpy2: %s \n time: %1.2e"%(resNumpy2,end)

print "-"*30
start = T.time()
resNumpy3 = combiProdNumpy3(a,b)
end = T.time() - start
print "Numpy2: %s \n time: %1.2e"%(resNumpy3,end)
