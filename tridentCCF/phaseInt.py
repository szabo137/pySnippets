"""
module to compute the general phase integral in CCF with airy functions
"""
import numpy as np
from scipy.special import airy
import matplotlib.pylab as plt

def b(e,c1,c2):
    return 2.0*np.pi*np.exp(-1j*e)/((3.0*c2)**(1.0/3.0))

def eta(r,c1,c2):
    return -r*c1/3.0/c2 + 2.0*c1**3/(27*c2**2)

def mu(r,c1,c2):
    return (r-c1**2/(3.0*c2))/((3.0*c2)**(1.0/3.0))


def B0ccf(r,c1=0.0,c2=0.0):
    et = eta(r,c1,c2)
    prefac = b(et,c1,c2)
    arg = mu(r,c1,c2)
    return prefac*airy(arg)[0]

def B1full(r,c1,c2):
    pass

if __name__=='__main__':
    args = np.linspace(0.0,20,300)
    print args.shape
    func=lambda x: B0ccf(x,1.5,1.3)*B0ccf(3-x,0.5,0.3)
    integFunc = lambda x: (func(x) - func(-x))/x
    #vals = np.array(map(func,args))
    vals=func(args)
    print vals.shape

    plt.plot(args,vals)
    plt.show()
