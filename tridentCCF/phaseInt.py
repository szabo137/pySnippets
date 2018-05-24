"""
module to compute the general phase integral in CCF with airy functions
"""
import numpy as np
from scipy.special import airy
import matplotlib.pylab as plt
from cQuad import cQuad
a=200


def b(e,c1,c2):
    return 2.0*np.pi*np.exp(1j*e)/((3.0*c2)**(1.0/3.0))

def eta(r,c1,c2):
    return -r*c1/3.0/c2 + 2.0*c1**3/(27*c2**2)

def mu(r,c1,c2):
    return (r-c1**2/(3.0*c2))/((3.0*c2)**(1.0/3.0))


def B0ccf(r,c1=0.0,c2=0.0):
    et = eta(r,c1,c2)
    prefac = b(et,c1,c2)
    arg = mu(r,c1,c2)
    return prefac*airy(arg)[0]

def B1ccf(r,c1,c2):
    et = eta(r,c1,c2)
    prefac = b(et,c1,c2)
    arg = mu(r,c1,c2)
    fac1 = c1/3.0/c2*airy(arg)[0]
    fac2 =1j/((3.0*c2)**(1.0/3.0))*airy(arg)[1] 
    return -prefac*(fac1 + fac2)

def B2ccf(r,c1,c2):
    et = eta(r,c1,c2)
    prefac = b(et,c1,c2)
    arg = mu(r,c1,c2)
    fac1 = ((c1/3.0/c2)**2 - arg/((3.0*c2)**(2.0/3.0)))*airy(arg)[0]
    fac2 =1j*2.0*c1/((3.0*c2)**(4.0/3.0))*airy(arg)[1] 
    return prefac*(fac1 + fac2)

def B0reg(r,c1,c2):
    return -1.0/r*(2.0*c1*B1ccf(r,c1,c2) + 3.0*c2*B2ccf(r,c1,c2))

def f(rStar):
    #return (rStar +15.2666666667
    #)/4.0970436
    return rStar


def integrandB1ccf(x,c1,c2):
    return x*np.exp(1j*c1*x**2/2.0 + 1j*c2*x**3/3.0)

def integB1ccf(r,c1,c2):
    func = lambda x: integrandB1ccf(x,c1,c2)
    args = np.linspace(-a,a,10000)
    vals = func(args)
    res = np.fft.fft(vals)
    
    return res


if __name__=='__main__':
    
    
    c1=-2.92894124925
    c2=0.976313749749
    
    args = np.linspace(1e-8,30,2000)
    print args.shape
    func0=lambda x: B0ccf(f(x),0.976313749749,0.325437916583)#*B1ccf(6.89928383156-f(x),-2.92894124925,0.976313749749)
    #func0=lambda x: B0ccf(x,c1,c2)
    func0reg=lambda x: B0reg(f(x),0.976313749749,0.325437916583)#*B1ccf(6.89928383156-f(x),-2.92894124925,0.976313749749)
    #func0reg=lambda x: B0reg(x,c1,c2)
    integFunc0 = lambda x: np.real((func0(x) - func0(-x))/x)
    integFunc0reg = lambda x: np.real((func0reg(x) - func0reg(-x))/x)
    #vals = np.array(map(func,args))
    #vals0=integFunc0(args)
    vals0=func0(args)
    print"test: %s"%(func0(0))
    #vals0reg=integFunc0reg(args)
    vals0reg=func0reg(args)
    #print vals.shape
    print vals0reg/vals0
    plt.plot(args,vals0,label='analytic')
    plt.plot(args,vals0reg,'-.',label='reg.')
    
    funcR=lambda x: 1/x
    #plt.plot(args,funcR(args),'.')
    funcB=lambda x: (2.0*c1*B1ccf(x,c1,c2) + 3.0*c2*B2ccf(x,c1,c2))
    #plt.plot(args,funcB(args),'x')
    print funcB(0)
    plt.legend()
    plt.show()
    """
    #fft check
    c1=-2.92894124925
    c2=0.976313749749
    #print B1ccf(1.3,c1,c2)
    args = np.linspace(-a,a,10000)
    integB1ccf(0.0,c1,c2)
    """
