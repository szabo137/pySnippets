"""
test calculation of momenta
"""


import numpy as np
import qft
import itertools
import time as T

def mom(E1,cTh1):
    rho=np.sqrt(E1**2-np.ones(len(E1)))
    return qft.MinkowskiVector([E1,rho*np.sqrt(np.ones(len(cTh1))-cTh1**2),0,rho*cTh1])


E1 = np.linspace(1.1,10.0,100)
cTh1 = np.linspace(0.0,0.9,100)


start = T.time()
mksMom = mom(E1,cTh1)
end = T.time() - start

print "mksMom (single): %s \n(time: %1.2e)"%(mksMom,end)


def mom2(e1,c1,e2,c2):
    E1,E2,cTh1,cTh2 = np.meshgrid(e1,e2,c1,c2)
    rho1=np.sqrt(E1**2-np.ones(E1.shape))
    rho2=np.sqrt(E2**2-np.ones(E2.shape))
    return np.array([qft.MinkowskiVector([E1,rho1*np.sqrt(np.ones(cTh1.shape)-cTh1**2),0,rho1*cTh1]),qft.MinkowskiVector([E2,rho2*np.sqrt(np.ones(cTh2.shape)-cTh2**2),0,rho2*cTh2])])


E1 = np.linspace(1.1,10.0,10)
cTh1 = np.linspace(0.0,0.9,10)
E2 = np.linspace(1.1,10.0,10)
cTh2 = np.linspace(0.0,0.9,10)



start = T.time()
res = mom2(E1,cTh1,E2,cTh2)
end = T.time() - start
print "mksMom (mesh): %s \n(time: %1.2e count:%s)"%(res,end,len(E1)*len(cTh1)*len(E2)*len(cTh2))
"""
def mom3(inArgs):
    e1 = inArgs[0]
    c1 = inArgs[1]
    e2 = inArgs[2]
    c2 = inArgs[3]
    rho1=np.sqrt(e1**2-1.0)
    rho2=np.sqrt(e2**2-1.0)
    return np.array([qft.MinkowskiVector([e1,rho1*np.sqrt(1.0-c1**2),0,rho1*c1]),qft.MinkowskiVector([e2,rho2*np.sqrt(1.0-c2**2),0,rho2*c2])])

start = T.time()
args=[t for t in itertools.product(E1,cTh1,E2,cTh2)]
#args=itertools.product(E1,cTh1,E2,cTh2)
res = map(mom3,args)
end = T.time() - start
print "mksMom (map): %s \n(time: %1.2e count:%s)"%(res,end,len(E1)*len(cTh1)*len(E2)*len(cTh2))
"""
