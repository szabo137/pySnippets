"""
routine to complex quadrature using scipy quad
"""
import scipy
from scipy.integrate import quad

def cQuad(func, a, b, **kwargs):
    def real_func(x,*args,**kwargs):
        return scipy.real(func(x,*args,**kwargs))
    def imag_func(x,*args,**kwargs):
        return scipy.imag(func(x,*args,**kwargs))
    print kwargs
    real_integral = quad(real_func, a, b, **kwargs)
    imag_integral = quad(imag_func, a, b, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0], real_integral[1:], imag_integral[1:])

if __name__=='__main__':
    import numpy as np
    testFunc = lambda x,a,b: a*x**2 + b + 1j*x**2
    res = cQuad(testFunc,-1,1,args=(4,1))
    print res
