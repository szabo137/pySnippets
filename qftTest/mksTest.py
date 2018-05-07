"""
test for mks - Vectors
"""
import numpy as np
import qft


E=np.linspace(1,3,3)
print "E: %s"%E
cTh = np.linspace(0.0,1.0,5)
print "cTh: %s"%cTh

def buildGrid(en,cth):
    rho=np.sqrt(en**2 - np.ones(en.shape))
    sinComp = np.sqrt(np.ones(cth.shape)-cth**2)
    return qft.MinkowskiVector([en,sinComp*rho,np.zeros(en.shape),cth*rho])


print"test %s"%buildGrid(np.array([1,2]),np.array([0.25,0.5]))


EE ,CC = np.meshgrid(E,cTh)

mom1= buildGrid(EE,CC)

EE2 ,CC2 = np.meshgrid(2*E,cTh*0.5)

mom2= buildGrid(EE2,CC2)

print mom1*mom2
