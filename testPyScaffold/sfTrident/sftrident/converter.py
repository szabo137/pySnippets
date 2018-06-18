"""
converts constants
"""
import math as m

def alpha2charge(alpha):
    r"""Converts finestructure constant to charge

    Parameters
    ----------
    alpha : float
        value of finestructure constant

    Returns
    -------
    float
        value of charge
    """
    return m.sqrt(4.0*m.pi*alpha)

def charge2alpha(charge):
    r"""Converts charge to finestructure constant

    Parameters
    ----------
    charge : float
        value of charge

    Returns
    -------
    float
        value of finestructure constant
    """
    return charge**2/(4.0*m.pi)
