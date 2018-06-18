import sys,os
sys.path.append(os.pardir+os.sep)

from Konstanten import *
from Grid_cms import *
from NME import *
from numpy import *


def eps_laser(l):
    return (l==1)*MinkowskiVector( [0,1,0,0] ) + (l==2)*MinkowskiVector( [0,0,1,0] )

def eps_photon(l):
    return (l==1)*MinkowskiVector([0,cos(theta_ph)*cos(phi_ph),cos(theta_ph)*sin(phi_ph),-sin(theta_ph)])\
           +(l==2)*MinkowskiVector( [0,-sin(phi_ph),cos(phi_ph),0])


w_laser     = 0.75/m
w_photon    = w_laser
w_cms       = w_laser
w           = w_laser
sqrt_s      = 2.*w_laser

pplusp = w_laser

a0            = 1.0
sigma         = 25 #pulse length
w_laser_l     = 1.e0
w_photon_l    = 4.e0
ksi           = pi/4.

Genau_theta   = 5
Genau_gamma   = 5

Harmonic_max  = 1


Envelope      = 'cos^2'

m_eff         = m*sqrt(1+0.5*a0**2)


theta_ph      = 0.
phi_ph        = 0.
phi_pos       = 0.


Abstand   = 2.0

y0_min  = lambda l: log ( w_cms/m - sqrt( w_cms**2/m**2 - 1./l ) )
y1_max  = lambda l: log ( w_cms/m + sqrt( w_cms**2/m**2 - 1./l ) )


p_perp=parray(1.0)
#print "pT: %s"%(p_perp)
#m_perp  = sqrt(1 + p_perp**2)
#p_perp  = parray(Grid_rap(y_z[0,:], w_cms, 0 , Harmonic_max + Abstand , Genau_gamma , a0)) #p_T
m_perp  = sqrt(1 + p_perp**2)


#y_z    = parray(linspace(y0_min(Harmonic_max + Abstand),y1_max(Harmonic_max + Abstand),Genau_theta),axis=1)[newaxis,:]
y_z = parray(np.log(pplusp/m_perp))

#print "yz: %s"%(y_z)
