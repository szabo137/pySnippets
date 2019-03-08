"""
transverse to spherical coordinates
"""
import numpy as np

def E(y,pt,m=1):
    mt = np.sqrt(pt**2 + m**2)
    #print "+++ %s"%mt
    return mt*np.cosh(y)

def c(y,pt,m=1):
    mt = np.sqrt(pt**2 + m**2)
    rho = np.sqrt(E(y,pt,m)**2 - m**2)
    #print "rho: %s"%rho
    #print "sh(y): %s"%np.sinh(y)
    return np.sinh(y)*mt/rho

def y(E,c,m=1):
    rho = np.sqrt(E**2 - m**2)
    pl = rho*c
    return .5*np.log((E+pl)/(E-pl))

def pt(E,c,m=1):
    rho = np.sqrt(E**2 - m**2)
    #if c>1.0:
        #print "## %s"%E
        #print "## %s"%c
    return rho*np.sin(np.arccos(c))

if __name__=='__main__':

    def trans2sphAndBack(Y,PT):
        tempE = E(Y,PT)
        tempC = c(Y,PT)
        #print "-"*20
        #print tempE
        #print tempC
        #print Y
        #print PT
        tempY = y(tempE,tempC)
        tempPT = pt(tempE,tempC)
        return tempY,tempPT

    testY = np.linspace(-4.0,4.0,200)
    testPT = np.linspace(1e-4,3.0,200)
    eps=1e-6
    #print c(0.6315789473684212,1e-7)
    #print np.sqrt(1e-8**2 +1)
    #print np.arccos(1.0000000000000002)
    for elY in testY:
        for elPT in testPT:
            yb,ptb = trans2sphAndBack(elY,elPT)
            if (np.abs(elY-yb)>eps) or (np.abs(elPT-ptb)>eps):
                print "y:   %s"%elY
                print "pt:  %s"%elPT
                print "yb:  %s"%yb
                print "ptb: %s"%ptb

    #print c(1.4,2.033333)
