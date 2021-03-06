# -*- coding: utf-8 -*-
"""
library for light cone coordinate formulation of momenta

Notes
-----
the kinparas are: ss,p2x,p2y,p2minus,p3x,p3y,p3minus

.. todo::
    -   implement preFac in getKin
    -   generalize the particle momenta, e.g. general mom(px,py,pminus)
"""
import numpy as np
import qft as qft


def laserOmega(ss):
    """
    lab energy of background photon from cms energy

    Parameters
    ----------
    ss : array_like
        center-of-momentum energy of the inital particles (one photon)

    Returns
    -------
    array_like
        energy of the background photon in the electron restsystem

    """
    return (ss**2 - 1)/(2.0)

def momPlus(momX,momY,momMinus):
    r"""plus component of onshell momentum

    Parameters
    ----------
    momX : float
        x component of the momentum

    momY : float
        y component of the momentum

    momMinus : float
        minus component of the momentum

    Returns
    -------
    float
        plus component of the momentum

    Notes
    -----
    Uses the onshell relation with unit mass:

    .. math:: p_+ = \frac{p_x^2 + p_y^2 + 1}{4p_-}
    """
    return (momX**2 + momY**2 + 1)/(4.0*momMinus)

def momLaser(ss):
    r"""
    momentum initial laser photon k

    Parameters
    ----------
    ss : float
        center-of-momentum energy of the inital particles (single photon)

    Returns
    -------
    qft.MinkowskiVector
        momentum of the inital background photon (electron rest system)

    Notes
    -----
    The momentum belongs to the z axis:

    .. math:: \left(k^\mu\right) = \left(w(s),0,0,w(s)\right)

    """
    omega = laserOmega(ss)
    return qft.MinkowskiVector([omega,0,0,omega])

def momInitElek():
    r"""
    momentum initial electron p

    Parameters
    ----------
    None

    Returns
    -------
    qft.MinkowskiVector
        momentum of inital electron (rest system)

    Notes
    -----
    Uses unity mass: :math:`\left(p^\mu \right) = (1,0,0,0)`
    """
    return qft.MinkowskiVector([1,0,0,0])

def momFinalElec1(p2x,p2y,p2minus):
    r"""
    momentum final electron no. 1: p_2

    Parameters
    ----------
    p2x : float
        x component of the momentum
    p2y : float
        y component of the momentum
    p2minus : float
        minus component of the momentum

    Returns
    -------
    qft.MinkowskiVector
        momentum four vector in cartesian coordinates

    Notes
    -----
        -   calculates the plus component with <momPlus>, e.g. uses onshell relation with unit mass
        -   calculates the 0th and 3rd component from light cone coordinates

            .. math:: p_0 = p_+ + p_-

            .. math:: p_3 = p_+ - p_-

    """
    p2plus = momPlus(p2x,p2y,p2minus)
    return qft.MinkowskiVector([p2plus + p2minus,p2x,p2y,p2plus - p2minus])

def momFinalPos(p3x,p3y,p3minus):
    r"""
    momentum final positron: p_3

    Parameters
    ----------
    p3x : float
        x component of the momentum
    p3y : float
        y component of the momentum
    p3minus : float
        minus component of the momentum

    Returns
    -------
    qft.MinkowskiVector
        momentum four vector in cartesian coordinates

    Notes
    -----
        -   calculates the plus component with <momPlus>, e.g. uses onshell relation with unit mass
        -   calculates the 0th and 3rd component from light cone coordinates

            .. math:: p_0 = p_+ + p_-

            .. math:: p_3 = p_+ - p_-

    """
    p3plus = momPlus(p3x,p3y,p3minus)
    return qft.MinkowskiVector([p3plus + p3minus,p3x,p3y,p3plus - p3minus])

def momFinalElec2(p2x,p2y,p2minus,p3x,p3y,p3minus):
    r"""
    momentum final electron no. 2: p_1

    Parameters
    ----------
    p2x : float
        x component of the momentum p2
    p2y : float
        y component of the momentum p2
    p2minus : float
        minus component of the momentum p2
    p3x : float
        x component of the momentum p3
    p3y : float
        y component of the momentum p3
    p3minus : float
        minus component of the momentum p3

    Returns
    -------
    qft.MinkowskiVector
        momentum four vector in cartesian coordinates

    Notes
    -----
        -   calculates the light cone components from energy momentum conservation:

        .. math:: p_{1x} = -p_{2x} - p_{3x}
        .. math:: p_{1y} = -p_{2y} - p_{3y}
        .. math:: p_{1-} = 0.5 -p_{2-} - p_{3-}

        -   calculates the plus component with <momPlus>, e.g. uses onshell relation with unit mass
        -   calculates the 0th and 3rd component from light cone coordinates

            .. math:: p_0 = p_+ + p_-

            .. math:: p_3 = p_+ - p_-

    """
    p1x = -p2x - p3x
    p1y = -p2y - p3y
    p1minus = 0.5-p2minus - p3minus
    p1plus = momPlus(p1x,p1y,p1minus)
    return qft.MinkowskiVector([p1plus + p1minus,p1x,p1y,p1plus - p1minus])


def charFuncPS(p1minus,p2minus,p3minus):
    r"""Characteristic function for minus components

    Parameters
    ----------
    p1minus : float
        minus component of the momentum p1

    p2minus : float
        minus component of the momentum p2

    p3minus : float
        minus component of the momentum p3

    Returns
    -------
    boolian
        only True if all input parameters are positiv

    """
    return (p1minus>0.0) and (p2minus>0.0) and (p3minus>0.0)

def preFacMatONSHELL(Plaser,Pelec,Pelec2):
    r"""Prefactor of the onshell matrix element

    Parameters
    ----------
    Plaser : qft.MinkowskiVector
        momentum of the laser photon

    Pelec : qft.MinkowskiVector
        momentum of the final electron no. 1

    Pelec2 : qft.MinkowskiVector
        momentum of the final electron no. 2

    Returns
    -------
    float
        prefactor of the onshell matrix element

    Notes
    -----
        -   direct and exchange part are selected within the input, e.g. 'exchange' means exchange of Pelec and Pelec2.
        -   calculates

        .. math:: -\frac{-i\pi}{k(p_{e1} - p_{e2})}

    """
    #only onshell!
    return -np.pi*1j/(2.0*(Plaser*(Pelec - Pelec2)))

def preFacMatOFFSHELL(Plaser,PinitElec,PfinalElec):
    r"""Prefactor of the offshell matrix element

    Parameters
    ----------
    Plaser : qft.MinkowskiVector
        momentum of the laser photon

    PinitElec : qft.MinkowskiVector
        momentum of the inital electron

    PfinalElec : qft.MinkowskiVector
        momentum of one final electron

    Returns
    -------
    float
        prefactor of the offshell matrix element

    Notes
    -----
        -   direct and exchange part are selected within the final electron, e.g. 'dirext' means inputting p1 and 'exchange' means inputting p2.
        -   calculates

        .. math:: \frac{1}{2k(p - p_{f})}

    """
    return 1.0/2.0/((PinitElec-PfinalElec)*Plaser)


def preFacRate(ss,p1m,p2m,p3m):
    r"""Prefactor of the full rate

    Parameters
    ----------

    ss : float
        center-of-momentum energy :math:`\sqrt(s)`

    p1m : float
        minus component of the momentum p1

    p2m : float
        minus component of the momentum p2

    p3m : float
        minus component of the momentum p3

    Returns
    -------
    float
        prefactor of the full rate

    Notes
    -----
        -   while this is of the level of the rate: there is no exchange part
        -   calculates

        .. math:: \frac{1}{4\omega^2}\frac{1}{(2\pi)^3 p_{2-}}\frac{1}{(2\pi)^3 p_{3-}}\frac{1}{p_{1-}}
    """
    omegaLaser = laserOmega(ss)
    return 1.0/(4.0*omegaLaser**2)/((2.0*np.pi)**3*2.0*p2m)/((2.0*np.pi)**3*2.0*p3m)/p1m

def getKin(kinPara):
    r"""Main function: Calculates the whole kinematic for a given parameter set

    Parameters
    ----------
    kinPara : list
        list of kinematic parameters [ss,p2x,p2y,p2minus,p3x,p3y,p3minus]

    Returns
    -------
    list
        list of [physArea,momenta,preFacs]

    Raises
    ------
    ValueError
        if the input is not a seven  parameter list

    Notes
    -----
        -   kinematic parameters are:
                -   'ss' center-of-momentum energy,
                -   'p2x'/'p2y'/'p2minus' light cone coordinates of p2,
                -   'p3x'/'p3y'/'p3minus' light cone coordinates of p3
        -   'physArea': boolian
        -   'momenta': [k,-p,p1,p2,p3]
        -   'prefacs': ['onshell direct','onshell exchange','offshell direct','offshell exchange','rate']
        -   if 'physArea' are False, return is [False,None,None]
    """
    try:
        [ss,p2x,p2y,p2minus,p3x,p3y,p3minus] = kinPara
    except ValueError as vE:
        raise ValueError("Not the correct kinPara format! <%s> given.\n(%s)"%(kinPara,vE))


    p1minus = 0.5-p2minus - p3minus # think about avoiding this double eval!
    physArea = charFuncPS(p1minus,p2minus,p3minus)
    if physArea:
        K=momLaser(ss)
        P=-momInitElek() #neg sign for direct usage in kinClass -> make this better?!
        P1 = momFinalElec2(p2x,p2y,p2minus,p3x,p3y,p3minus)
        P2 = momFinalElec1(p2x,p2y,p2minus)
        P3 = momFinalPos(p3x,p3y,p3minus)
        preFac1 = preFacMatONSHELL(K,-P,P1) # direct onshell
        preFac2 = preFacMatONSHELL(K,-P,P2) # exchange onshell
        preFac3 = preFacMatOFFSHELL(K,-P,P1) # direct offshell
        preFac4 = preFacMatOFFSHELL(K,-P,P2) # exchange offshell
        preFac5 = preFacRate(ss,p1minus,p2minus,p3minus)
        return [physArea,np.array([K,P,P1,P2,P3]),[preFac1,preFac2,preFac3,preFac4,preFac5]]
    else:
        return [physArea,None,None]
