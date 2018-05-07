"""
enviroment to test spinor products with gamma matrices
"""
import numpy as np
import qft
import time

gamma=qft.GammaMatrix()
gammaVec = np.array([gamma[item].data for item in np.arange(4)])
g5=qft.gamma5()


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

def buildMomSingle(E,cTh):
    rho = np.sqrt(E**2 - 1.0)
    return qft.MinkowskiVector([E,rho*np.sqrt(1.0-cTh*cTh),0,cTh*rho])



def prodFOR(spinorA,spinorB):
    #current1 = np.zeros(4,dtype=np.complex)
    #current2 = np.zeros(4,dtype=np.complex)
    #current3 = np.zeros(4,dtype=np.complex)
    #current4 = np.zeros(4,dtype=np.complex)
    current1 = np.zeros(4,dtype=object)
    current2 = np.zeros(4,dtype=object)
    current3 = np.zeros(4,dtype=object)
    current4 = np.zeros(4,dtype=object)
    factor=g5
    for item in np.arange(4):
        diracMatrix = gamma[item]*factor
        current1[item] = (spinorA[0]*(diracMatrix*spinorB[0]))
        current2[item] = (spinorA[1]*(diracMatrix*spinorB[1]))
        current3[item] = (spinorA[0]*(diracMatrix*spinorB[1]))
        current4[item] = (spinorA[1]*(diracMatrix*spinorB[0]))
    return np.array([qft.MinkowskiVector(current1),qft.MinkowskiVector(current2),qft.MinkowskiVector(current3),qft.MinkowskiVector(current4)])


def prodTOB(spinorA,spinorB):
    factor = g5
    return np.array(map(qft.MinkowskiVector,np.swapaxes(np.array(spinorA[np.newaxis,np.newaxis,:]*gamma[np.newaxis,:,np.newaxis]*factor*spinorB[:,np.newaxis,np.newaxis]),1,2).reshape(4,4)))



if __name__=='__main__':
    energy = np.random.uniform(1.0,10.0,10)#np.linspace(1.0,10.0,10)
    cosTheta = np.random.uniform(-1.0,1.0,10)#np.linspace(-1.0,1.0,10)
    
    #P=buildMomSingle(energy,cosTheta)
    P=buildMom(energy,cosTheta)
    
    UBars = np.array([qft.SpinorUBar((P,1.0),spin) for spin in [1,2]])
    Vs = np.array([qft.SpinorV((P,1.0),spin) for spin in [1,2]])
    
    print "mom: %s"%P
    print "ubar: \n%s"%UBars
    print "v: \n%s"%Vs
    
    #gammaVec = np.array([gamma[item].data for item in np.arange(4)])
    #print gammaVec.shape
    
    start = time.time()
    resFOR = prodFOR(UBars,Vs)
    end = time.time() -start
    print "resFOR: \n%s"%(resFOR)
    print "time: %1.2e"%end

    start = time.time()
    resTOB = prodTOB(UBars,Vs)
    end = time.time() -start
    print "resTOB: \n%s"%(resTOB)
    """
    for item,el in enumerate(resTOB):
        for item2,el2 in enumerate(el):
            for item3,el3 in enumerate(el2):
                print "(%s,%s,%s): %s"%(item,item2,item3,el3)
    """
    print "time: %1.2e"%end
