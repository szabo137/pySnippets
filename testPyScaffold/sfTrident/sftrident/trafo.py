"""
Transformation from spherical and transverse coordinates to lightcone
"""
import numpy as np
import numdifftools as nd

def spherical2lightcone(E,cTh,phi,m=1.0):
    rho = np.sqrt(E**2 - m**2)*(E>=1)
    sTh = np.sqrt(1-cTh**2)
    return ((E-rho*cTh)/2.0,rho*sTh*np.cos(phi),rho*sTh*np.sin(phi))

def transverse2lightcone(y,pT,phi,m=1.0):
    return (0.5*np.exp(-y)*np.sqrt(pT**2 + m**2),pT*np.cos(phi),pT*np.sin(phi))

def trafoSingle(ss,E2,cTh2,phi2,y3,pT3,phi3,m2=1.0,m3=1.0):
    p2minus,p2x,p2y = spherical2lightcone(E2,cTh2,phi2,m2)
    p3minus,p3x,p3y = transverse2lightcone(y3,pT3,phi3,m3)
    #print p2minus*np.sqrt(E2**2-1)*p3minus*pT3
    return np.array([ss,p2x,p2y,p2minus,p3x,p3y,p3minus])

def trafo(newCoordniates):
    return trafoSingle(*newCoordniates)

def jacobian(newCoordniates):
    ss,E2,cTh2,phi2,y3,pT3,phi3 = newCoordniates
    p2minus = spherical2lightcone(E2,cTh2,phi2)[0]
    p3minus = transverse2lightcone(y3,pT3,phi3)[0]
    return p2minus*np.sqrt(E2**2-1)*p3minus*pT3


def jacobianNUM(newCoordniates):
    jac = nd.Jacobian(trafo)
    return np.linalg.det(jac(newCoordniates))


if __name__=='__main__':
    import time
    ss=np.array([3.353])
    E=np.array([1.2])
    cTh = np.array([0.93])
    phi1=np.array([0.0])
    y = np.array([1.8])
    pT=np.array([0.4])
    phi2 = np.array([0.0])
    print trafo([ss,E,cTh,phi1,y,pT,phi2])
    start = time.time()
    jacNum = jacobianNUM([ss,E,cTh,phi1,y,pT,phi2])
    end = time.time() - start
    print "jacobian numerical: %s (time: %1.2e)"%(jacNum,end)

    start = time.time()
    jacAna = jacobian([ss,E,cTh,phi1,y,pT,phi2])
    end = time.time() - start
    print "jacobian analytic: %s (time: %1.2e)"%(jacAna,end)
