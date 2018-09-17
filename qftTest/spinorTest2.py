"""
tests the relations of spinors in qft lib
"""
import numpy as np
import qft

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




if __name__=='__main__':
    energy = np.array([3.5])#np.random.uniform(1.0,10.0,10)#np.linspace(1.0,10.0,10)
    cosTheta = np.array([0.9])#np.random.uniform(-1.0,1.0,10)#np.linspace(-1.0,1.0,10)

    P=buildMom(energy,cosTheta)
    print "mom: %s"%(P)
    print "mom: %s"%(str((P.data).T))
    print "mass: %s"%(P*P)


    S = qft.SpinorU((P,1.0),1)
    SBar = qft.SpinorUBar((P,1.0),1)
    print "U(p)"
    print "sp: %s"%(S)
    print "spBar: %s"%(SBar)
    print "spBar*sp: %s"%(SBar*S)

    V = qft.SpinorV((P,1.0),1)
    VBar = qft.SpinorVBar((P,1.0),1)
    print "V(p)"
    print "sp: %s"%(V)
    print "spBar: %s"%(VBar)
    print "spBar*sp: %s"%(VBar*V)




    V2 = qft.SpinorU((-P,1.0),1)
    VBar2 = qft.SpinorUBar((-P,1.0),1)
    print "U(-p)"
    print "sp: %s"%(V2)
    print "spBar: %s"%(VBar2)
    print "spBar*sp: %s"%(VBar*V)


    energy = np.array([2.8])#np.random.uniform(1.0,10.0,10)#np.linspace(1.0,10.0,10)
    cosTheta = np.array([0.8])#np.random.uniform(-1.0,1.0,10)#np.linspace(-1.0,1.0,10)
    P2=buildMom(energy,cosTheta)
    print "mom: %s"%(P)
    print "mom: %s"%(str((P.data).T))
    print "mass: %s"%(P*P)



    gamma = qft.GammaMatrix()

    print qft.SpinorUBar((P,1.0),1)*(gamma[0]*qft.SpinorU((P2,1.0),1))
    print qft.SpinorUBar((P,1.0),1)*(gamma[0]*qft.SpinorV((-P2,1.0),1))
