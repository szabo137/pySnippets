"""
test enviroment to integrate over a delta function
"""
import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import quad


def lorentzCurve(x,eps):
    return eps/(x**2 + eps**2)/np.pi

testEps = 1e-13
x0=3.0
integFunc = lambda t: lorentzCurve(t-x0,testEps)*(t+1)
res = quad(integFunc,x0-300*testEps,x0+300*testEps,points=([x0]))

args = np.linspace(-10,10,300)
vars = integFunc(args)
plt.plot(args,vars)
plt.title("$\\epsilon=%s$\n$integ=%s$"%(testEps,res))
plt.show()
