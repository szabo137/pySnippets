"""
kinematic utility functions
"""
"""
contains the polarisation vectors of photon lines
"""

import qft

def BGpolarisationBase(n):
    if n==1:
        return qft.MinkowskiVector([0.0,1.0,0.0,0.0])
    elif n==2:
        return qft.MinkowskiVector([0.0,0.0,1.0,0.0])
    else:
        raise ValueError("Index of polarisation base vector needs to be 1 or 2 (<%s> given)"%n)



def photoNumBW(r,mom):
    """
    returns the general photoNum for the Breit Wheeler parts: s_r

    does not change for the exchange part

    todo:
        -   insert in kinClass and return as Function
    """
    #for item,el in enumerate(mom):
    #    print "mom %s: %s"%(item,el)
    pt = mom[2] + mom[3] + mom[4]
    return (pt*pt - 1)/(2.0*(-mom[0]*mom[1])) - r

def photoNumConshell(mom):
    """
    returns the onshell photo number for the compton vertex: r_on

    momFormat: k,-p,p1
    """
    #for item,el in enumerate(mom):
    #    print"internal photoNum mom %s: %s"%(item,el)
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
    return (l==1)*qft.MinkowskiVector( [0,1,0,0] ) + (l==2)*qft.MinkowskiVector( [0,0,1,0] )
