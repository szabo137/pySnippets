"""
module to test the evaluation of airy functions
"""
from scipy.special import airy
import numpy as np
import time as T
import matplotlib.pylab as plt


def AiryEval(z):
    return airy(z)[0]

rounds = 100
args=np.linspace(-10,5,rounds)

start = T.time()
vals=AiryEval(args)
end = T.time() - start

plt.plot(args,vals)
plt.title("Airy function - full time: %1.2e s - avg. time: %1.2e s (%s rounds)"%(end,end/rounds,rounds))

plt.show()
