"""
lib for light cone coordinate formulation of momenta

Notes
-----
the kinparas are: ss,p2x,p2y,p2minus,p3x,p3y,p3minus

TODO:
    -   implement preFac in getKin
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
    """plus component of onshell momentum

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
    Uses the onshell relation with unit mass

    .. math:: p_\mu p^\mu = 1

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
    Uses unity mass:

    .. math:: \left(p^\mu \right) = (1,0,0,0)
    
    """
    return qft.MinkowskiVector([1,0,0,0])

def momFinalElec1(p2x,p2y,p2minus):
    """
    momentum final electron no. 1: p_2

    Parameters
    ----------
    p2x : float
        x component of the momentum

    """
    p2plus = momPlus(p2x,p2y,p2minus)
    return qft.MinkowskiVector([p2plus + p2minus,p2x,p2y,p2plus - p2minus])

def momFinalPos(p3x,p3y,p3minus):
    """
    momentum final positron p_3
    """
    p3plus = momPlus(p3x,p3y,p3minus)
    return qft.MinkowskiVector([p3plus + p3minus,p3x,p3y,p3plus - p3minus])

def momFinalElec2(p2x,p2y,p2minus,p3x,p3y,p3minus):
    p1x = -p2x - p3x
    p1y = -p2y - p3y
    p1minus = 0.5-p2minus - p3minus
    p1plus = momPlus(p1x,p1y,p1minus)
    return qft.MinkowskiVector([p1plus + p1minus,p1x,p1y,p1plus - p1minus])


def charFuncPS(p1minus,p2minus,p3minus):
    return (p1minus>0.0) and (p2minus>0.0) and (p3minus>0.0)

def preFacMatONSHELL(Plaser,Pelec,Pelec2):
    #only onshell!
    return -np.pi*1j/(2.0*(Plaser*(Pelec - Pelec2)))

def preFacMatOFFSHELL(Plaser,PinitElec,PfinalElec):
    return 1.0/2.0/((PinitElec-PfinalElec)*Plaser)


def preFacRate(ss,p1m,p2m,p3m):
    #full prefac
    omegaLaser = laserOmega(ss)
    return 1.0/(4.0*omegaLaser**2)/((2.0*np.pi)**3*2.0*p2m)/((2.0*np.pi)**3*2.0*p3m)/p1m

def getKin(kinPara):
    """
    returns the right format for kinClass and proofs physArea
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
