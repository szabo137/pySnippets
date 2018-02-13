"""
sum the elements of an n x m array to one n-array
"""

import numpy as np
from itertools import izip, imap
import time
#testArr=np.array([[0,0,0,0],[0,1,1,0],[0,1,2,0],[0,0,1,0]])
testArr = np.arange(10000).reshape(100,100)
start = time.time()
#rowSum=[sum(el) for el in testArr]
#rowSum=np.array(map(sum,testArr))
#rowSum=np.array([t for t in imap(np.sum,testArr)])
rowSum=np.sum(testArr,axis=1)
end = time.time() -start

print"row sum: %s (time: %1.2e)"%(rowSum,end)
print type(rowSum)

start = time.time()
#colSum = map(np.sum,zip(*testArr))
colSum = np.array(map(np.sum,testArr.T))
#colSum = np.array([t for t in imap(sum,izip(*testArr))])
end = time.time() -start

print "column sum: %s (time: %1.2e)"%(colSum,end)
print type(colSum)
