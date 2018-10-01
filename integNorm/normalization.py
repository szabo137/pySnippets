"""
module to compute the normalisation of a 2D contour plot
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec
from itertools import imap
from scipy.integrate import dblquad

# Simple data to display in various forms
def fkt(a,b):
    return (a**2 + b**2)*(np.abs(b)<2)*(np.abs(a)<2)

resDirect = dblquad(fkt,-1, 1, lambda x: -1, lambda x: 1) #=8/3
print resDirect

def flattenContour(Z,axis='row',func=np.sum):
    if axis=='row':
        return np.array(map(func,Z))
    elif axis=='col':
        return np.array(map(func,Z.T))
    else:
        raise ValueError("axis need to be 'row' or 'col'. ('%s' given)"%axis)


x = np.linspace(0,2, 10)
y = np.linspace(0,1, 10)

XX,YY = np.meshgrid(x,y)
ZZ=fkt(XX,YY)

print ZZ
print XX[XX<1]

CS = plt.contourf(XX,YY,ZZ)
cbar=plt.colorbar(CS,format="%1.2e")
plt.show()
