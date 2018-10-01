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
    
    print "---------- XFEL set -----------"
    Eelec = 50.0/restMass
    photoEn = omegaFromE(Eelec,cmsEn)
    print "sqrt(s) = %s me (%s MeV)"%(cmsEn,cmsEn*restMass)
    print "init Eelec = %s me (%s MeV)"%(Eelec,Eelec*restMass)
    print "init omega = %s me (%s MeV)"%(photoEn,photoEn*restMass)
    print "relative rapidity = %s "%(rapidity(photoEn,Eelec))

    print "---------- e^- rest -----------"
    Eelec = 1.0
    photoEn = omegaFromE(Eelec,cmsEn)
    print "sqrt(s) = %s me (%s MeV)"%(cmsEn,cmsEn*restMass)
    print "init Eelec = %s me (%s MeV)"%(Eelec,Eelec*restMass)
    print "init omega = %s me (%s MeV)"%(photoEn,photoEn*restMass)
    print "relative rapidity = %s "%(rapidity(photoEn,Eelec))

    print "---------- equal mom. -----------"
    Eelec = (cmsEn**2 + 1)/(2*cmsEn)
    photoEn = omegaFromE(Eelec,cmsEn)
    print "sqrt(s) = %s me (%s MeV)"%(cmsEn,cmsEn*restMass)
    print "init Eelec = %s me (%s MeV)"%(Eelec,Eelec*restMass)
    print "init omega = %s me (%s MeV)"%(photoEn,photoEn*restMass)
    print "relative rapidity = %s "%(rapidity(photoEn,Eelec))

    print "---------- LUXE -----------"
    Eelec = 17.5*1e3/restMass
    photoEn = omegaFromE(Eelec,cmsEn)
    print "sqrt(s) = %s me (%s MeV)"%(cmsEn,cmsEn*restMass)
    print "init Eelec = %s me (%s MeV)"%(Eelec,Eelec*restMass)
    print "init omega = %s me (%s MeV)"%(photoEn,photoEn*restMass)
    print "relative rapidity = %s "%(rapidity(photoEn,Eelec))

    print "---------- E-144 -----------"
    Eelec = 46.6*1e3/restMass
    photoEn = omegaFromE(Eelec,cmsEn)
    print "sqrt(s) = %s me (%s MeV)"%(cmsEn,cmsEn*restMass)
    print "init Eelec = %s me (%s MeV)"%(Eelec,Eelec*restMass)
    print "init omega = %s me (%s MeV)"%(photoEn,photoEn*restMass)
    print "relative rapidity = %s "%(rapidity(photoEn,Eelec))
