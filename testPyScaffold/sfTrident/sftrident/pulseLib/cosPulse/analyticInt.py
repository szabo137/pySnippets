"""
contains the analytic expressions of the internal integrals for cos pulses

used the bounds: (0,phi) from Seipt
"""
import numpy as np
from math import pi
#from sfTrident.settings import laserConstants
#dphi = laserConstants['dphi']/np.pi


def __I1(phi,dphi):
    """
    analytic solution of

        integal(cos^2(phi/(2dphi))*cos(phi), phi)

    from wolfram alpha
    """
    dphi=dphi/pi
    term1 = np.sin(phi)*(dphi**2*np.cos(phi/dphi) + dphi**2 - 1.0)
    term2 = -dphi*np.cos(phi)*np.sin(phi/dphi)
    return (term1 + term2)/(2.0*(dphi**2 - 1.0))


def __I2(phi,dphi):
    """
    analytic solution of

        integal(cos^2(phi/(2dphi))*sin(phi), phi)

    from wolfram alpha
    """
    dphi=dphi/pi
    term1 = np.cos(phi)*(dphi**2*np.cos(phi/dphi) + dphi**2 - 1.0)
    term2 = dphi*np.sin(phi)*np.sin(phi/dphi)
    return -(term1 + term2)/(2.0*(dphi**2 - 1.0))


def __I3(phi,dphi):
    """
    analytic solution of

        integal(cos^4(phi/(2dphi))*cos^2(phi), phi)

    from wolfram alpha
    """
    dphi=dphi/pi
    term1 = 8.0*dphi/(2.0*dphi+1.0)*np.sin((1.0/dphi+2.0)*phi)
    term2 = 16.0*dphi*np.sin(phi/dphi)
    term3 = 2.0*dphi*np.sin(2.0*phi/dphi)
    term4 = dphi/(dphi - 1.0)*np.sin(2*(dphi-1.0)/dphi*phi)
    term5 = dphi/(dphi+1.0)*np.sin(2.0*(dphi+1.0)/dphi*phi)
    term6 = 8.0*dphi/(2.0*dphi - 1.0)*np.sin((2*dphi-1.0)/dphi*phi)
    term7 = 12.0*phi
    term8 = 6.0*np.sin(2.0*phi)
    return (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8)/64.0

def __I4(phi,dphi):
    """
    analytic solution of

        integal(cos^4(phi/(2dphi))*sin^2(phi), phi)

    from wolfram alpha
    """
    dphi=dphi/pi
    term1 = -8.0*dphi/(2.0*dphi+1.0)*np.sin((1.0/dphi+2.0)*phi)
    term2 = 16.0*dphi*np.sin(phi/dphi)
    term3 = 2.0*dphi*np.sin(2.0*phi/dphi)
    term4 = -dphi/(dphi - 1.0)*np.sin(2*(dphi-1.0)/dphi*phi)
    term5 = -dphi/(dphi+1.0)*np.sin(2.0*(dphi+1.0)/dphi*phi)
    term6 = -8.0*dphi/(2.0*dphi - 1.0)*np.sin((2*dphi-1.0)/dphi*phi)
    term7 = 12.0*phi
    term8 = -6.0*np.sin(2.0*phi)
    return (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8)/64.0

def analyticInternalInt(item):
    """
        returns the internal integral I_item as a function of phi
        (the returned functions need to have an array as first arg)

        todo: - generalze the bounds
    """
    integrals = [lambda x,dphi: __I1(x,dphi) - __I1(np.zeros(len(x)),dphi),
                lambda x,dphi: __I2(x,dphi) - __I2(np.zeros(len(x)),dphi),
                lambda x,dphi: __I3(x,dphi) - __I3(np.zeros(len(x)),dphi),
                lambda x,dphi: __I4(x,dphi) - __I4(np.zeros(len(x)),dphi)]
    return integrals[item]
    #return lambda x: integrals[item](dphi) - integrals[item](x)
