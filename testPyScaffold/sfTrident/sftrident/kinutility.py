# -*- coding: utf-8 -*-
r"""
kinematic utility functions
"""

import qft

def BGpolarisationBase(n):
    r"""Polarisation base of the background field

    Parameters
    ----------
    n : int
        index of element

    Returns
    -------
    qft.MinkowskiVector
        polariation base vector of the given index

    Raises
    ------
    ValueError
        if <n> is not 1 or 2

    Notes
    -----
    The polarisation base vectors are:

    .. math:: \epsilon^\mu_1 = (0,1,0,0)

    and

    .. math:: \epsilon^\mu_2 = (0,0,1,0)

    .. todo::
        -   remove this function and use <laserPolarisation>
    """
    if n==1:
        return qft.MinkowskiVector([0.0,1.0,0.0,0.0])
    elif n==2:
        return qft.MinkowskiVector([0.0,0.0,1.0,0.0])
    else:
        raise ValueError("Index of polarisation base vector needs to be 1 or 2 (<%s> given)"%n)



def photoNumBW(r,mom):
    r"""
    general photoNum for the Breit-Wheeler parts: math:'s_r'

    Parameters
    ----------
    r : float
        photo number of the compton part
    mom : np.ndarray(shape=(5,))
        set of all momenta (k,-p,p1,p2,p3)

    Returns
    -------
    float
        photo number of the Breit-Wheeler part (general form)

    Notes
    -----
    The formular for the general Breit-Wheeler photo number is

    .. math:: s_r = \frac{p_t^2 - 1}{2(kp)} - r

    where :math:`p_t = p_1 + p_2 + p_3` and the mass is set to unity.
        -   This does not change for the exchange part!
        -   The

    .. todo::
        -   insert in kinClass and return as Function
    """
    pt = mom[2] + mom[3] + mom[4]
    return (pt*pt - 1)/(2.0*(-mom[0]*mom[1])) - r

def photoNumConshell(mom):
    r"""
    photoNum for the Compton onshell part: math:'r_{\mathrm{on}}'

    Parameters
    ----------
    mom : np.ndarray(shape=(3,))
        set of compton momenta: (k,-p,p1) or (k,-p,p2)

    Returns
    -------
    float
        photo number of the Compton onshell part

    Notes
    -----
    The formular for the general Breit-Wheeler photo number is

    .. math:: r_{\mathrm{on}} = -\frac{\delta p}{2(\delta p \cdot k)} - r

    where :math:`\delta p = p-p1` or :math:`\delta p = p-p2`, respectivly and the mass is set to unity.

    .. todo::
        -   insert in kinClass and return as Function
    """
    deltaP = -mom[2]-mom[1]
    return -(deltaP*deltaP)/(2.0*deltaP*mom[0])

def rFromRstar(rStar,mom):
    """
    trafo r -> rStar

    use CMomenta
    """
    deltaP = -mom[2]-mom[1]
    return (rStar-(deltaP*deltaP))/(2.0*deltaP*mom[0])

def laserPolarisation(l):
    r"""Polarisation base of the background field

    Parameters
    ----------
    l : int
        index of element

    Returns
    -------
    qft.MinkowskiVector
        polariation base vector of the given index

    Raises
    ------
    ValueError
        if <l> is not 1 or 2

    Notes
    -----
    The polarisation base vectors are:

    .. math:: \epsilon^\mu_1 = (0,1,0,0)

    and

    .. math:: \epsilon^\mu_2 = (0,0,1,0)

    """
    return (l==1)*qft.MinkowskiVector( [0,1,0,0] ) + (l==2)*qft.MinkowskiVector( [0,0,1,0] )
