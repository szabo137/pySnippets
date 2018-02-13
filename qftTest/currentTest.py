"""
tests structures to evaluate products like ubar*gamma*v etc
"""
import numpy as np
import qft
import time as T
gamma=qft.GammaMatrix()
def buildMom(E,cTh):
    if not(isinstance(E,np.ndarray) and isinstance(cTh,np.ndarray)):
        if np.isscalar(E):
            E=[E]
        if np.isscalar(cTh):
            cTh=[cTh]
        E=np.array(E)
        cTh = np.array(cTh)
    rho = np.sqrt(E**2 - np.ones(E.shape))
    return qft.MinkowskiVector([E,rho*np.sqrt(np.ones(cTh.shape)-cTh*cTh),0,cTh*rho])

def buildSpinorU(mom,m=1.0):
    return np.array([qft.SpinorU((mom,m),s) for s in [1,2]])

def buildSpinorUBar(mom,m=1.0):
    return np.array([qft.SpinorUBar((mom,m),s) for s in [1,2]])

def buildSpinorV(mom,m=1.0):
    return np.array([qft.SpinorV((mom,m),s) for s in [1,2]])

def buildSpinorVBar(mom,m=1.0):
    return np.array([qft.SpinorVBar((mom,m),s) for s in [1,2]])

def buildProd(spinorA, spinorB):
    current1 = np.zeros(4,dtype=np.complex)
    current2 = np.zeros(4,dtype=np.complex)
    for item in np.arange(4):
        #print "full %s"%(complex(spinorA[0]*(gamma[item]*spinorB[0])))
        #print "type %s"%(type(spinorA[0]*(gamma[item]*spinorB[0])))
        current1[item] = (spinorA[0]*(gamma[item]*spinorB[0]))[0]
        current2[item] = (spinorA[1]*(gamma[item]*spinorB[1]))[0]
    print type(current1)
    return np.array([qft.MinkowskiVector(current1),qft.MinkowskiVector(current2)])
    



e=1.2
c=0.8

P=buildMom(e,c)
print "P: %s"%P
print "P type: %s"%(type(P))
print "P shape: %s"%(str(P._0().shape))
U=buildSpinorU(P)
print "U: %s"%(str(U.shape))
print "U[newaxis,:]: %s"%(str(U[None,:].shape))
print "U[newaxis,newaxis,:]: %s"%(str(U[None,None,:].shape))
print "U: %s"%(str(U.shape))
print "U[:,newaxis]: %s"%(str(U[:,None].shape))
print "U[:,newaxis,newaxis]: %s"%(str(U[:,None,None].shape))
print "U[newaxis,:,newaxis]: %s"%(str(U[None,:,None]))
Ubar=buildSpinorUBar(P)


start = T.time()
res = buildProd(Ubar,U)
end= T.time() - start
#res = (gammaVec[newaxis,:,newaxis]*Ubar[:,newaxis,newaxis])
print "res: %s"%(res)
print "type: %s"%(type(res))
print "shape: %s"%(str(res.shape))
print "time: %1.2e"%(end)

for el in res:
    print "test: %s"%(el*P)
    print "test type: %s"%(type(el*P))
    print "test pass?: %s"%(el*P==2.0)
