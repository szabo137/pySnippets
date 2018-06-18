"""
current functions from tobias (copy, selfimplemet)
"""
import numpy as np
import sftrident.qft as qft

m=1.0
a0=0.5
xi=np.pi/4.0
theta_ph      = 0.
phi_ph        = 0.
def eps_photon(l):
    return (l==1)*qft.MinkowskiVector([0,np.cos(theta_ph)*np.cos(phi_ph),np.cos(theta_ph)*np.sin(phi_ph),-np.sin(theta_ph)])\
           +(l==2)*qft.MinkowskiVector( [0,-np.sin(phi_ph),np.cos(phi_ph),0])

def eps_laser(l):
    return (l==1)*qft.MinkowskiVector( [0,1,0,0] ) + (l==2)*qft.MinkowskiVector( [0,0,1,0] )

def J00(upb,vpp,epsp):
    res = upb*(qft.feyndagg(epsp)*vpp)
    return res

def J1plus(upb,vpp,epsp,k,pe,pp):
    #pol = epsPM(2)
    pol= eps_laser(1)*np.cos(xi)-1j*eps_laser(2)*np.sin(xi)
    dp=m*a0/(4.0*(k*pe))
    dpp=m*a0/(4.0*(k*pp))
    term1 = upb*qft.feyndagg(pol)*qft.feyndagg(k)*qft.feyndagg(epsp)*vpp
    term2 = upb*qft.feyndagg(epsp)*qft.feyndagg(k)*qft.feyndagg(pol)*vpp
    return dp*term1 - dpp*term2 # eigl -

def J1minus(upb,vpp,epsp,k,pe,pp):
    pol = eps_laser(1)*np.cos(xi)+1j*eps_laser(2)*np.sin(xi)
    dp=m*a0/(4.0*(k*pe))
    dpp=m*a0/(4.0*(k*pp))
    term1 = upb*(qft.feyndagg(pol)*(qft.feyndagg(k)*(qft.feyndagg(epsp)*vpp)))
    term2 = upb*qft.feyndagg(epsp)*qft.feyndagg(k)*qft.feyndagg(pol)*vpp
    return dp*term1 - dpp*term2 #mglw mit +

def J02(upb,vpp,epsp,k,pe,pp):
    dp=m*a0/(4.0*(k*pe))
    dpp=m*a0/(4.0*(k*pp))
    return -4.0*dp*dpp*(upb*qft.feyndagg(k)*vpp)*(k*epsp)


def J22(upb,vpp,epsp,k,pe,pp):
    dp=m*a0/(4.0*(k*pe))
    dpp=m*a0/(4.0*(k*pp))
    return -2.0*((np.cos(xi))**2-(np.sin(xi))**2)*dp*dpp*(upb*qft.feyndagg(k)*vpp)*(k*epsp)

def run():
    pPos = qft.MinkowskiVector([ 1.41518852,  1. ,         0. ,         0.05252185])
    pEl = qft.MinkowskiVector([ 1.41518852 ,-1.   ,       0.      ,    0.05252185])
    pKlaser = qft.MinkowskiVector([ 1.46771037 , 0.   ,       0.    ,     -1.46771037])
    pKphoton = qft.MinkowskiVector([ 1.46771037 , 0.    ,      0.    ,      1.46771037])
    #polMinus = qft.MinkowskiVector([ 0.00000000+0.j     ,     0.70710678+0.j       ,   0.00000000-0.70710678j,  0.00000000+0.j        ])
    #polPlus = qft.MinkowskiVector([ 0.00000000+0.j      ,    0.70710678+0.j         , 0.00000000+0.70710678j,  0.00000000+0.j        ])



    V     = np.array([qft.SpinorV ((pPos, m),s1) for s1 in [1,2]])
    Ubar  = np.array([qft.SpinorUBar((pEl,m),s1) for s1 in [1,2]])
    print "photoPol1: %s"%(eps_photon(l=1))
    print "photoPol2: %s"%(eps_photon(l=2))
    print "curr pol=0: %s"%(str(J1plus(Ubar[0],V[0],eps_photon(l=1),pKlaser,pEl,pPos)))
    print "curr pol=1: %s"%(str(J1plus(Ubar[0],V[0],eps_photon(l=2),pKlaser,pEl,pPos)))

if __name__=='__main__':
    pPos = qft.MinkowskiVector([ 1.41518852,  1. ,         0. ,         0.05252185])
    pEl = qft.MinkowskiVector([ 1.41518852 ,-1.   ,       0.      ,    0.05252185])
    pKlaser = qft.MinkowskiVector([ 1.46771037 , 0.   ,       0.    ,     -1.46771037])
    pKphoton = qft.MinkowskiVector([ 1.46771037 , 0.    ,      0.    ,      1.46771037])
    #polMinus = qft.MinkowskiVector([ 0.00000000+0.j     ,     0.70710678+0.j       ,   0.00000000-0.70710678j,  0.00000000+0.j        ])
    #polPlus = qft.MinkowskiVector([ 0.00000000+0.j      ,    0.70710678+0.j         , 0.00000000+0.70710678j,  0.00000000+0.j        ])



    V     = np.array([qft.SpinorV ((pPos, m),s1) for s1 in [1,2]])
    Ubar  = np.array([qft.SpinorUBar((pEl,m),s1) for s1 in [1,2]])
    print J1plus(Ubar[0],V[0],eps_photon(l=1),pKlaser,pEl,pPos)
