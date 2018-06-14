"""
internal integrals in constant crossed field approx -> requires xi==0.0
"""
import numpy as np

def __I1(phi,dphi):
    """
        int f_1(phi) dphi
        
    with f_1(phi) = phi
    """
    return phi**2/2.0

def __I3(phi,dphi):
    """
        int f_3(phi) dphi
        
    with f_3(phi) = f_1(phi)^2
    """
    return phi**3/3.0

def __I2(phi,dphi):
    """
    just a token, is not needed in ccf
    """
    return 0.0

def __I4(phi,dphi):
    """
    just a token, is not needed in ccf
    """
    return 0.0

def internalInt(item):
    """
        returns the internal integral I_item as a function of phi
        (the returned functions need to have an array as first arg)

        todo: - generalze the bounds
    """
    integrals = [lambda x,dphi: __I1(x,dphi),
                lambda x,dphi: __I2(x,dphi),
                lambda x,dphi: __I3(x,dphi) ,
                lambda x,dphi: __I4(x,dphi)]
    return integrals[item]
    #return lambda x: integrals[item](dphi) - integrals[item](x)
