"""
calculation of the rapidity shift for heads on collisions
"""
import numpy as np


def omegaFromE(E,ss):
    """
    calculation of omega from given init electron energy and cms energy
    """
    rho = np.sqrt(E**2-1)
    return (ss**2 - 1)/(2*(E+rho))


def rapidity(omega,E):
    """
    calculation of relative rapidity from lab to moving frame
    """
    Erel = omega + E
    PLrel = omega - np.sqrt(E**2 -1)
    return 0.5*np.log((Erel + PLrel)/(Erel - PLrel))


if __name__=='__main__':
    restMass = 0.511 #MeV
    cmsEn = 3.353
    Eelec = 50.0/restMass
    photoEn = omegaFromE(Eelec,cmsEn)
    print "sqrt(s) = %s MeV"%(cmsEn*restMass)
    print "init Eelec = %s MeV"%(Eelec*restMass)
    print "init omega = %s MeV"%(photoEn*restMass)
    print "relative rapidity = %s "%(rapidity(photoEn,Eelec))
