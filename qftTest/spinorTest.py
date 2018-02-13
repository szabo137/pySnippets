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
    energy = np.random.uniform(1.0,10.0,10)#np.linspace(1.0,10.0,10)
    cosTheta = np.random.uniform(-1.0,1.0,10)#np.linspace(-1.0,1.0,10)
    
    P=buildMom(energy,cosTheta)
    print "mom: %s"%(P)
    print "mom: %s"%(str((P.data).T))
    print "mass: %s"%(P*P)
    
    test = map(qft.MinkowskiVector,(P.data).T)
    print test
    
    S = qft.SpinorU((P,1.0),1)
    SBar = qft.SpinorUBar((P,1.0),1)
    print "sp: %s"%(S)
    print "spBar: %s"%(SBar)
    print "spBar*sp: %s"%(SBar*S)

    V = qft.SpinorV((P,1.0),1)
    VBar = qft.SpinorVBar((P,1.0),1)
    print "sp: %s"%(V)
    print "spBar: %s"%(VBar)
    print "spBar*sp: %s"%(VBar*V)

    Us = np.array([qft.SpinorU((P,1.0),1),qft.SpinorU((P,1.0),2)])
    #print "U1, U2: %s"%(Us)
    
    UBars = np.array([qft.SpinorUBar((P,1.0),1),qft.SpinorUBar((P,1.0),2)])
    #print "U1bar, U2bar: %s"%(UBars)
    
    res= np.sum(Us*UBars)
    gamma = qft.GammaMatrix()
    expectResU = qft.feyndagg(P)+1.0
    
    print "res-expect: %s"%(np.abs((res - expectResU).data).all()<=1e-10)    
    Vs = np.array([qft.SpinorV((P,1.0),1),qft.SpinorV((P,1.0),2)])
    #print "V1, V2: %s"%(Vs)
    
    VBars = np.array([qft.SpinorVBar((P,1.0),1),qft.SpinorVBar((P,1.0),2)])
    #print "V1bar, V2bar: %s"%(VBars)
    
    res= np.sum(Vs*VBars)
    gamma = qft.GammaMatrix()
    expectResV = qft.feyndagg(P)-1.0
    
    S = qft.SpinorU((P,1.0),1)
    res3=expectResV*S
    print "EOM U: %s"%(np.abs(res3.data).all()<=1e-10)
    
    S = qft.SpinorV((P,1.0),1)
    res3=expectResU*S
    print "EOM U: %s"%(np.abs(res3.data).all()<=1e-10)
    
    S = qft.SpinorUBar((P,1.0),1)
    res3=S*expectResV
    print "EOM U: %s"%(np.abs(res3.data).all()<=1e-10)
    
    S = qft.SpinorVBar((P,1.0),1)
    res3=S*expectResU
    print "EOM U: %s"%(np.abs(res3.data).all()<=1e-10)
    
