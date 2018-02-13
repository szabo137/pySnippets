"""
test the products of spinors and gamma matrices
"""
import numpy as np
import qft
import time as T

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

start = T.time()
E,cTh = np.linspace(1.2,3.0,10),np.linspace(-1.0,1.0,10)
EE,CC = np.meshgrid(E,cTh)
P1 = buildMom(EE,CC)
end=T.time() - start
print "build mom: %1.2e s"%end
P2 = buildMom(1.3,0.6)
#print 2*P1
start = T.time()
U1bar = qft.SpinorUBar((P1,1.0),1)
U2=qft.SpinorU((P1,1.0),1)
end=T.time() - start
print "build spinor: %1.2e s"%end
gamma=qft.GammaMatrix()
#comp0 = np.complex(U1bar*(gamma[0]*U2))
#comp1 = np.complex(U1bar*(gamma[1]*U2))
#comp2 = np.complex(U1bar*(gamma[2]*U2))
#comp3 = np.complex(U1bar*(gamma[3]*U2))
start = T.time()
comp0 = U1bar*(gamma[0]*U2)
comp1 = U1bar*(gamma[1]*U2)
comp2 = U1bar*(gamma[2]*U2)
comp3 = U1bar*(gamma[3]*U2)
end=T.time() - start
print "build comp: %1.2e s"%end
start = T.time()
res = qft.MinkowskiVector([comp0,comp1,comp2,comp3])
end=T.time() - start
print "build mks: %1.2e s"%end

start = T.time()
rest = res*P1
end=T.time() - start
print "build prod: %1.2e s"%end
