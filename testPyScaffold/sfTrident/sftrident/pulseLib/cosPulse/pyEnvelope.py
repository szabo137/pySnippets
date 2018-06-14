"""
contains the cos square envelope
"""
import numpy as np
#from sfTrident.settings import laserConstants
#dphi = laserConstants['dphi']

def envelope(phi,dphi):
    return (np.cos(np.pi*phi/(2.0*dphi)))**2*(phi<= dphi)*(phi>=-dphi)
