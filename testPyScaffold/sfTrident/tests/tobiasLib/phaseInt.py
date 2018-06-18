"""
return phaseIntegrals from tobias lib
"""
import numpy as np
import NME as qft

from Config_paar import *
from A_Funktion import *
from Kinematik import *
import numpy
from zerhacker import *
from NME import *
from asimps import *

m=1.0

def fd(x):
    return feyndagg(x)

def tobPI():
    p_pos,p_el,k_laser,k_photon,q_pos,eps_m,eps_p = kinematik()
    momenta = [p_pos,p_el,k_laser,k_photon,q_pos,eps_m,eps_p]
    para = [a0,ksi,sigma]
    if Envelope == 'cos^2':
        a = sigma/w_laser
    elif Envelope == 'Gauss':
        a = sigma*5.
    elif Envelope == 'Box':
        a = sigma+1


    #@zerhacker( slicearguments = (0,1) , slicelength = 250 , axis = 1)
    #@zerhacker( slicearguments = (0,1) , slicelength = 250 , axis = 0)
    def integrate(p_el,p_pos,k_photon,k_laser,M,N):
        return asimps( lambda x: A_m_n_nSVEA(M,N,x,p_el,p_pos,k_photon,k_laser),- a, a, Nx=421 , errorabs = 1e-4, maxrecur = 100 )


    A_11  = integrate(p_el,p_pos,k_photon,k_laser,1,1)
    #print "A11: %s"%(str(A_11))
    A_1_1 = integrate(p_el,p_pos,k_photon,k_laser,1,-1)
    #print "A1-1: %s"%(str(A_1_1))
    A_20  = integrate(p_el,p_pos,k_photon,k_laser,2,0)
    #print "A20: %s"%(str(A_20))


    A_22  = integrate(p_el,p_pos,k_photon,k_laser,2,2)
    #print "A22: %s"%(str(A_22))
    A_2_2  = integrate(p_el,p_pos,k_photon,k_laser,2,-2)

    A_00  = A_0_0(A_11[0],A_1_1[0],A_20[0],A_22[0],A_2_2[0] )
    #print "A00: %s"%(str(A_00))

    return (A_11[0],A_1_1[0],A_20[0],A_22[0],A_2_2[0],A_00,momenta,para)


def phaseIntFromTob():


    #calc tobias phase ints
    A11,A1_1,A20,A22,A2_2,A00,momenta,para = tobPI()
    [pPos,pEl,pKlaser,pKphoton,q_pos,eps_m,eps_p] = momenta
    xi=para[1]
    #calc my phase ints from tobias
    B1 = np.cos(xi)/2.0*pKlaser.minus()/2.0*(A11 + A1_1)
    B2 = np.sin(xi)/(2.0*1j)*pKlaser.minus()/2.0*(A11 - A1_1)
    B3 = pKlaser.minus()/4.0*A20 + pKlaser.minus()/8.0*np.cos(2.0*xi)*(A22 + A2_2)
    B0 = pKlaser.minus()/2.0*A00
    return (B1,B2,B3,B0,momenta,para)
